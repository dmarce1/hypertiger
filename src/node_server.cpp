/*
 * node_server.cpp
 *
 *  Created on: Jun 11, 2015
 *      Author: dmarce1
 */

#include "defs.hpp"
#include "node_server.hpp"
#include "future.hpp"
#include "options.hpp"

#include <array>
#include <streambuf>
#include <fstream>
#include <iostream>
#include <sys/stat.h>
#if !defined(_MSC_VER)
#include <unistd.h>
#endif

#include <hpx/include/lcos.hpp>
#include <hpx/include/util.hpp>
#include <hpx/runtime/threads/run_as_os_thread.hpp>

std::vector<real> node_server::outflows(physics::NF, 0.0);

real node_server::get_rotation_count() const {
	return current_time;
}

hpx::future<void> node_server::exchange_flux_corrections() {
	const geo::octant ci = my_location.get_child_index();
	constexpr auto full_set = geo::face::full_set();
	for (auto& f : full_set) {
		const auto face_dim = f.get_dimension();
		auto const& this_aunt = aunts[f];
		if (!this_aunt.empty()) {
			std::array<integer, NDIM> lb, ub;
			lb[XDIM] = lb[YDIM] = lb[ZDIM] = 0;
			ub[XDIM] = ub[YDIM] = ub[ZDIM] = INX;
			if (f.get_side() == geo::MINUS) {
				lb[face_dim] = 0;
			} else {
				lb[face_dim] = INX;
			}
			ub[face_dim] = lb[face_dim] + 1;
			auto data = grid_ptr->get_flux_restrict(lb, ub, face_dim);
			this_aunt.send_hydro_flux_correct(std::move(data), f.flip(), ci);
		}
	}

	constexpr integer size = geo::face::count() * geo::quadrant::count();
	std::array<hpx::future<void>, size> futs;
	for (auto& f : futs) {
		f = hpx::make_ready_future();
	}
	integer index = 0;
	for (auto const& f : geo::face::full_set()) {
		if (this->nieces[f]) {
			for (auto const& quadrant : geo::quadrant::full_set()) {
				futs[index++] =
						niece_hydro_channels[f][quadrant].get_future().then(
								hpx::util::annotated_function(
										[this, f, quadrant](hpx::future<std::vector<real> > && fdata) -> void
										{
											const auto face_dim = f.get_dimension();
											std::array<integer, NDIM> lb, ub;
											switch (face_dim) {
												case XDIM:
												lb[XDIM] = f.get_side() == geo::MINUS ? 0 : INX;
												lb[YDIM] = quadrant.get_side(0) * (INX / 2);
												lb[ZDIM] = quadrant.get_side(1) * (INX / 2);
												ub[XDIM] = lb[XDIM] + 1;
												ub[YDIM] = lb[YDIM] + (INX / 2);
												ub[ZDIM] = lb[ZDIM] + (INX / 2);
												break;
												case YDIM:
												lb[XDIM] = quadrant.get_side(0) * (INX / 2);
												lb[YDIM] = f.get_side() == geo::MINUS ? 0 : INX;
												lb[ZDIM] = quadrant.get_side(1) * (INX / 2);
												ub[XDIM] = lb[XDIM] + (INX / 2);
												ub[YDIM] = lb[YDIM] + 1;
												ub[ZDIM] = lb[ZDIM] + (INX / 2);
												break;
												case ZDIM:
												lb[XDIM] = quadrant.get_side(0) * (INX / 2);
												lb[YDIM] = quadrant.get_side(1) * (INX / 2);
												lb[ZDIM] = f.get_side() == geo::MINUS ? 0 : INX;
												ub[XDIM] = lb[XDIM] + (INX / 2);
												ub[YDIM] = lb[YDIM] + (INX / 2);
												ub[ZDIM] = lb[ZDIM] + 1;
												break;
											}
											grid_ptr->set_flux_restrict(fdata.get(), lb, ub, face_dim);
										},
										"node_server::exchange_flux_corrections::set_flux_restrict"));
			}
		}
	}
	return hpx::when_all(std::move(futs));
}

void node_server::all_hydro_bounds() {
	exchange_interlevel_hydro_data();
	collect_hydro_boundaries();
//	send_hydro_amr_boundaries();
	++hcycle;
}

void node_server::exchange_interlevel_hydro_data() {

	if (is_refined) {
		for (auto const& ci : geo::octant::full_set()) {
			auto data = child_hydro_channels[ci].get_future(hcycle).get();
			grid_ptr->set_restrict(data, ci);
		}
	}
	auto data = grid_ptr->get_restrict();
	integer ci = my_location.get_child_index();
	if (my_location.level() != 0) {
		parent.send_hydro_children(std::move(data), ci, hcycle);
	}
}
void node_server::collect_hydro_boundaries() {
	for (auto const& dir : geo::direction::full_set()) {
		if (!neighbors[dir].empty()) {
			auto bdata = grid_ptr->get_hydro_boundary(dir);
			neighbors[dir].send_hydro_boundary(std::move(bdata), dir.flip(),
					hcycle);
		}
	}

	std::array<hpx::future<void>, geo::direction::count()> results;
	integer index = 0;
	for (auto const& dir : geo::direction::full_set()) {
		if (!(neighbors[dir].empty() /*&& my_location.level() == 0*/)) {
			results[index++] =
					sibling_hydro_channels[dir].get_future(hcycle).then(
							hpx::util::annotated_function(
									[this](hpx::future<sibling_hydro_type> && f) -> void
									{
										auto&& tmp = f.get();
										grid_ptr->set_hydro_boundary(tmp.data, tmp.direction);
									},
									"node_server::collect_hydro_boundaries::set_hydro_boundary"));
		}
	}
	wait_all_and_propagate_exceptions(std::move(results));

	for (auto& face : geo::face::full_set()) {
		if (my_location.is_physical_boundary(face)) {
			grid_ptr->set_physical_boundaries(face, current_time);
		}
	}
}

void node_server::send_hydro_amr_boundaries() {

	if (is_refined) {
		constexpr auto full_set = geo::octant::full_set();
		for (auto& ci : full_set) {
			const auto& flags = amr_flags[ci];
			for (auto& dir : geo::direction::full_set()) {
				if (flags[dir]) {
					std::array<integer, NDIM> lb, ub;
					const integer width = BW;
					get_boundary_size(lb, ub, dir, OUTER, INX, width);
					for (integer dim = 0; dim != NDIM; ++dim) {
						lb[dim] = ((lb[dim] - BW)) + 2 * BW
								+ ci.get_side(dim) * (INX);
						ub[dim] = ((ub[dim] - BW)) + 2 * BW
								+ ci.get_side(dim) * (INX);
					}
					auto data = grid_ptr->get_prolong(lb, ub);
					children[ci].send_hydro_boundary(std::move(data), dir,
							hcycle);
				}
			}
		}
	}
}

inline bool file_exists(const std::string& name) {
	struct stat buffer;
	return (stat(name.c_str(), &buffer) == 0);
}

//HPX_PLAIN_ACTION(grid::set_omega, set_omega_action2);
//HPX_PLAIN_ACTION(grid::set_pivot, set_pivot_action2);

std::size_t node_server::load_me(FILE *fp, bool old_format) {
	std::size_t cnt = 0;
	auto foo = std::fread;
	refinement_flag = false;
	cnt += foo(&step_num, sizeof(integer), 1, fp) * sizeof(integer);
	cnt += foo(&current_time, sizeof(real), 1, fp) * sizeof(real);
	cnt += grid_ptr->load(fp);
	return cnt;
}

std::size_t node_server::save_me(std::ostream& strm) const {
	std::size_t cnt = 0;

	cnt += write(strm, step_num);
	cnt += write(strm, current_time);

	assert(grid_ptr != nullptr);
	cnt += grid_ptr->save(strm);
	return cnt;
}

void node_server::save_to_file(const std::string& fname,
		std::string const& data_dir) const {
	hpx::util::high_resolution_timer timer;
	save(0, data_dir + fname).get();
	double elapsed = timer.elapsed();
	printf("Saving took %f seconds\n", elapsed);
}

void node_server::load_from_file(const std::string& fname,
		std::string const& data_dir) {
	hpx::util::high_resolution_timer timer;

	integer rec_size = 0;
	int total_nodes;
	hpx::threads::run_as_os_thread([&]() {
		FILE* fp = fopen((data_dir + fname).c_str(), "rb");
		if (fp == NULL) {
			printf("Failed to open file\n");
			abort();
		}
		fseek(fp, -sizeof(integer), SEEK_END);
		std::size_t read_cnt = fread(&rec_size, sizeof(integer), 1, fp);
		fread(outflows.data(), sizeof(real), physics::NF, fp)*sizeof(real);
		fseek(fp, -4 * sizeof(real) - sizeof(integer), SEEK_END);
		fclose(fp);

		// work around limitation of ftell returning 32bit offset
			std::ifstream in((data_dir + fname).c_str(), std::ifstream::ate | std::ifstream::binary);
			std::size_t end_pos = in.tellg();

			total_nodes = end_pos / rec_size;
		}).get();

	printf("Loading %d nodes\n", total_nodes);

	load(0, total_nodes, rec_size, opts.output_only, data_dir + fname).get();
	double elapsed = timer.elapsed();
	printf("Loading took %f seconds\n", elapsed);
}

void node_server::load_from_file_and_output(const std::string& fname,
		const std::string& outname, std::string const& data_dir) {
	load_from_file(fname, data_dir);
}

void node_server::clear_family() {
	me = hpx::invalid_id;
	std::fill(aunts.begin(), aunts.end(), hpx::invalid_id);
	std::fill(nieces.begin(), nieces.end(), false);
}

integer child_index_to_quadrant_index(integer ci, integer dim) {
	integer index;
	if (dim == XDIM) {
		index = ci >> 1;
	} else if (dim == ZDIM) {
		index = ci & 0x3;
	} else {
		index = (ci & 1) | ((ci >> 1) & 0x2);
	}
	return index;
}

void node_server::initialize(real t) {
	hcycle = 0;
	step_num = 0;
	refinement_flag = 0;
	is_refined = false;
	neighbors.resize(geo::direction::count());
	nieces.resize(NFACE);
	aunts.resize(NFACE);
	current_time = t;
	dx = TWO * opts.xscale / real(INX << my_location.level());
	for (auto& d : geo::dimension::full_set()) {
		xmin[d] = opts.xscale * my_location.x_location(d);
	}
	if (current_time == ZERO) {
		grid_ptr = std::make_shared < grid > (dx, xmin, true);
	} else {
		grid_ptr = std::make_shared < grid > (dx, xmin);
	}
}

node_server::node_server(const node_location& loc, const node_client& parent_id,
		real t, std::size_t _step_num, std::size_t _hcycle) :
		my_location(loc), parent(parent_id) {
	initialize(t);
	step_num = _step_num;
	hcycle = _hcycle;
}

node_server::node_server(const node_location& _my_location, integer _step_num,
		bool _is_refined, real _current_time,
		const std::array<integer, NCHILD>& _child_d, grid _grid,
		const std::vector<hpx::id_type>& _c, std::size_t _hcycle) {
	my_location = _my_location;
	initialize(_current_time);
	hcycle = _hcycle;
	is_refined = _is_refined;
	step_num = _step_num;
	current_time = _current_time;
	grid_ptr = std::make_shared < grid > (std::move(_grid));
	if (is_refined) {
		std::copy(_c.begin(), _c.end(), children.begin());
	}
	child_descendant_count = _child_d;
}

void node_server::recv_hydro_boundary(std::vector<simd_vector>&& bdata,
		const geo::direction& dir, std::size_t cycle) {
	sibling_hydro_type tmp;
	tmp.data = std::move(bdata);
	tmp.direction = dir;
	sibling_hydro_channels[dir].set_value(std::move(tmp), cycle);
}

void node_server::recv_hydro_children(std::vector<real>&& data,
		const geo::octant& ci, std::size_t cycle) {
	child_hydro_channels[ci].set_value(std::move(data), cycle);
}

void node_server::recv_hydro_flux_correct(std::vector<real>&& data,
		const geo::face& face, const geo::octant& ci) {
	const geo::quadrant index(ci, face.get_dimension());
	niece_hydro_channels[face][index].set_value(std::move(data));
}

void node_server::amr_driver() {
	integer output_cnt;

	printf("Starting...\n");
	regrid(me.get_gid(), false);

	real output_dt = opts.output_dt;

	printf("output_dt = %e\n", output_dt);
	real& t = current_time;
	integer step_num = 0;

	auto fut_ptr = me.get_ptr();
	node_server* root_ptr = fut_ptr.get();

	output_cnt = root_ptr->get_rotation_count() / output_dt;

	real bench_start, bench_stop;

	while (current_time < opts.stop_time) {
		if (step_num > opts.stop_step)
			break;

		auto time_start = std::chrono::high_resolution_clock::now();
		if (!opts.disable_output
				&& root_ptr->get_rotation_count() / output_dt >= output_cnt) {
			//	if (step_num != 0) {

			std::string fname = "X." + std::to_string(int(output_cnt)) + ".chk";
			save_to_file(fname, opts.data_dir);
			printf("doing silo out...\n");

			fname = "X." + std::to_string(int(output_cnt));
			output(opts.data_dir, fname, output_cnt);
			++output_cnt;

		}
		if (step_num == 0) {
			bench_start = hpx::util::high_resolution_clock::now() / 1e9;
		}

		real dt = 0;

		integer next_step = (std::min)(step_num + 15, opts.stop_step + 1);
		real omega_dot = 0.0, omega = 0.0, theta = 0.0, theta_dot = 0.0;

		auto rc = step(next_step - step_num).get();
		for (integer f = 0; f != physics::NF; ++f) {
			outflows[f] += rc.first[f].sum();
		}
		dt = rc.second;

		auto sums = grid_ptr->sums();

		hpx::threads::run_as_os_thread([=]()
		{
			printf("!!! %i %e %e %e\n", int(next_step - 1), double(t), double(dt), sums[0].sum());
		});
		step_num = next_step;

		if (step_num % 15 == 0) {
			auto need_break = hpx::threads::run_as_os_thread([&]()
			{

				bench_stop = hpx::util::high_resolution_clock::now() / 1e9;
				return false;
			});
			if (need_break.get())
				break;
		}
	}
	bench_stop = hpx::util::high_resolution_clock::now() / 1e9;
	{
		if (!opts.disable_output)
			output(opts.data_dir, "final", output_cnt);
	}
}

hpx::future<std::vector<simd_vector>> node_server::grid_step() {
#if HPX_HAVE_ITTNOTIFY != 0 && !defined(HPX_HAVE_APEX)
	static hpx::util::itt::string_handle sh("node_server::nonrefined_step");
	hpx::util::itt::task t(hpx::get_thread_itt_domain(), sh);
#endif

	real cfl0 = cfl;
	dt_ = ZERO;

	all_hydro_bounds();

	grid_ptr->store();
	hpx::future<void> fut = hpx::make_ready_future();

	hpx::shared_future<real> dt_fut = global_timestep_channel.get_future();
	auto outflow = std::make_shared < std::vector
			< simd_vector >> (physics::NF, 0.0);

	for (integer rk = 0; rk < NRK; ++rk) {

		fut =
				fut.then(
						hpx::util::annotated_function(
								[rk, cfl0, this, dt_fut, outflow](hpx::future<void> f)
								{
									f.get();        // propagate exceptions

									real a = grid_ptr->compute_fluxes();
									hpx::future<void> fut_flux = is_refined ? hpx::make_ready_future() : exchange_flux_corrections();

									if (rk == 0) {
										const real dx = TWO * opts.xscale /
										real(INX << my_location.level());
										dt_ = cfl0 * dx / a;
										local_timestep_channels[NCHILD].set_value(dt_);
									}

									return fut_flux.then(
											hpx::launch::async(hpx::threads::thread_priority_boost),
											hpx::util::annotated_function(
													[rk, this, dt_fut, outflow](hpx::future<void> f)
													{
														f.get(); // propagate exceptions

														grid_ptr->compute_sources(current_time);
														grid_ptr->compute_dudt();

														if (rk == 0) {
															dt_ = dt_fut.get();
														}
														*outflow = grid_ptr->next_u(rk, current_time, dt_, std::move(*outflow));

														all_hydro_bounds();
														return outflow;

													}, "node_server::nonrefined_step::compute_fmm"
											));
								},
								"node_server::nonrefined_step::compute_fluxes"));
	}

	return fut.then(hpx::launch::sync, [this, outflow](hpx::future<void>&& f)
	{
		f.get(); // propagate exceptions...
			current_time += dt_;
			return *outflow;
		});
}

hpx::future<std::pair<std::vector<simd_vector>, real>> node_server::local_step(
		integer steps) {
	auto fut = hpx::make_ready_future(
			std::make_pair(std::vector<simd_vector>(physics::NF, 0.0), 0.0));
	for (integer i = 0; i != steps; ++i) {
		fut =
				fut.then(
						[this, i, steps](hpx::future<std::pair<std::vector<simd_vector>,real>> fut)
						{
							auto time_start = std::chrono::high_resolution_clock::now();
							auto next_dt = timestep_driver_descend();

							auto outflow = grid_step().get();

							auto tmp = fut.get();
							auto dt = tmp.second;
							for( integer f = 0; f != physics::NF; ++f) {
								outflow[f] += tmp.first[f];
							}
							if (my_location.level() == 0)
							{
								double time_elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(
										std::chrono::high_resolution_clock::now() - time_start).count();

								if (i + 1 != steps)
								{
									hpx::threads::run_as_os_thread([=]()
											{
												printf("%i %e %e %e\n", int(step_num), double(current_time), double(dt),
														time_elapsed);
											}); // do not wait for output to finish
								}
							}
							++step_num;

							return std::make_pair(outflow,next_dt.get());
						});
	}
	return fut;
}

hpx::future<std::pair<std::vector<simd_vector>, real>> node_server::step(
		integer steps) {
	grid_ptr->set_coordinates();

	std::array<hpx::future<void>, NCHILD> child_futs;
	if (is_refined) {
		for (integer ci = 0; ci != NCHILD; ++ci) {
			child_futs[ci] = children[ci].step(steps);
		}
	}

	auto fut = local_step(steps);

	if (is_refined) {
		return hpx::dataflow(hpx::launch::sync,
				[this](hpx::future<std::pair<std::vector<simd_vector>, real>> fut, hpx::future<void>&& f)
				{
					f.get(); // propagate exceptions
					return fut;
				}, std::move(fut), hpx::when_all(std::move(child_futs)));
	}

	return fut;
}

void node_server::timestep_driver_ascend(real dt) {
	global_timestep_channel.set_value(dt);
	if (is_refined) {
		for (auto& child : children) {
			child.timestep_driver_ascend(dt);
		}
	}
}
void node_server::set_local_timestep(integer idx, real dt) {
	local_timestep_channels[idx].set_value(dt);
}

hpx::future<real> node_server::timestep_driver_descend() {
	if (is_refined) {
		std::array<hpx::future<real>, NCHILD + 1> futs;
		integer index = 0;
		for (auto& local_timestep : local_timestep_channels) {
			futs[index++] = local_timestep.get_future();
		}

		return hpx::dataflow(hpx::launch::sync,
				hpx::util::annotated_function(
						[this](std::array<hpx::future<real>, NCHILD+1> dts_fut) -> double
						{
							auto dts = hpx::util::unwrapped(dts_fut);
							real dt = *std::min_element(dts.begin(), dts.end());

							if (my_location.level() == 0)
							{
								timestep_driver_ascend(dt);
							}
							else
							{
								parent.set_local_timestep(my_location.get_child_index(), dt);
							}

							return dt;
						}, "node_server::timestep_driver_descend"), futs);
	} else {
		return local_timestep_channels[NCHILD].get_future().then(
				hpx::launch::sync,
				[this](hpx::future<real>&& f)
				{
					real dt = f.get();
					parent.set_local_timestep(my_location.get_child_index(), dt);
					return dt;
				});
	}
}

/*
 * node_server_actions_1.cpp
 *
 *  Created on: Sep 23, 2016
 *      Author: dmarce1
 */

#include "defs.hpp"
#include "container_device.hpp"
#include "future.hpp"
#include "node_client.hpp"
#include "node_server.hpp"
#include "options.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <fstream>
#include <vector>

#include <hpx/include/lcos.hpp>
#include <hpx/include/run_as.hpp>
#include <hpx/lcos/broadcast.hpp>

#include <boost/iostreams/stream.hpp>

void parallel_output_gather(grid::output_list_type);
void parallel_output_complete(std::string dirname, std::string fname,
		int cycle);
HPX_PLAIN_ACTION(node_server::parallel_output_complete, parallel_output_complete_action);

std::stack<grid::output_list_type> node_server::pending_output;

void node_server::parallel_output_gather(grid::output_list_type&& list) {
	static hpx::mutex mtx;
	if (!list.nodes.empty()) {
		std::lock_guard<hpx::mutex> lock(mtx);
		pending_output.push(std::move(list));
	}
}

void node_server::parallel_output_complete(std::string dirname,
		std::string fname, real tm, int cycle) {
	grid::output_list_type olist;
	while (!pending_output.empty()) {
		auto next_list = std::move(pending_output.top());
		pending_output.pop();
		grid::merge_output_lists(olist, std::move(next_list));
	}
	grid::output(std::move(olist), dirname, fname, tm, cycle);

}

hpx::future<grid::output_list_type> node_server::load(integer cnt,
		integer total_nodes, integer rec_size, bool do_output,
		std::string filename) {
	if (my_location.level() == 0)
		me = hpx::invalid_id;
	else
		me = this->get_unmanaged_id();

	char flag = '0';
	std::vector<integer> counts(NCHILD);

	// run output on separate thread
	std::size_t read_cnt = 0;
	hpx::threads::run_as_os_thread([&]() {
		FILE* fp = fopen(filename.c_str(), "rb");
		fseek(fp, cnt * rec_size, SEEK_SET);
		read_cnt += fread(&flag, sizeof(char), 1, fp);

		for (auto& this_cnt : counts)
		{
			read_cnt += fread(&this_cnt, sizeof(integer), 1, fp);
		}
		//  printf( "RECSIZE=%i\n", int(rec_size));
			load_me(fp, rec_size==65739);
			fclose(fp);
		}).get();

	std::array<hpx::future<grid::output_list_type>, NCHILD> futs;
	if (flag == '1') {
		is_refined = true;

		integer index = 0;
		for (auto const& ci : geo::octant::full_set()) {
			integer loc_id = ((cnt * options::all_localities.size())
					/ (total_nodes + 1));

			futs[index++] =
					hpx::new_ < node_server
							> (options::all_localities[loc_id], my_location.get_child(
									ci), me.get_gid(), ZERO, step_num, hcycle).then(
									[this, ci, counts, do_output, total_nodes, rec_size, filename](hpx::future<hpx::id_type>&& fut)
									{
										children[ci] = fut.get();
										return children[ci].load(counts[ci], total_nodes, rec_size, do_output, filename);
									});

		}
	} else if (flag == '0') {
		is_refined = false;
		std::fill_n(children.begin(), NCHILD, node_client());
	} else {
		printf("Corrupt checkpoint file\n");
		//		sleep(10);
		hpx::this_thread::sleep_for(std::chrono::seconds(10));
		abort();
	}

	if (!is_refined && do_output) {
		auto l = grid_ptr->get_output_list();
		if (opts.parallel_silo) {
			parallel_output_gather(std::move(l));
		}
		return hpx::make_ready_future(l);
	}

	return hpx::dataflow(
			[this, do_output, filename](std::array<hpx::future<grid::output_list_type>, NCHILD>&& futs)
			{
				grid::output_list_type my_list;
				for (auto&& fut : futs) {
					if (fut.valid()) {
						if (do_output) {
							grid::merge_output_lists(my_list, fut.get());
						} else {
							fut.get();
						}
					}
				}
				if (my_location.level() == 0) {
					if (do_output) {
						auto silo_name = opts.output_filename;
						/* Skip for now, more interested in SILO */
						//	if (hydro_on && opts.problem == DWD) {
						//		diagnostics();
						//	}
						std::string this_fname;
						printf("Outputing...\n");
						if (opts.parallel_silo) {
							std::string dir_name = silo_name + std::string(".silo.data");
							if (system((std::string("mkdir -p ") + dir_name + std::string("\n")).c_str()) != 0) {
								abort_error();
							}
							const auto cycle = get_rotation_count();
							const auto sz = opts.all_localities.size();
							std::vector<hpx::future<void>> futs(sz);
							for (integer i = 0; i != sz; ++i) {
								this_fname = dir_name + std::string("/") + silo_name + std::string(".") + std::to_string(i) + std::string(".silo");
								futs[i] = hpx::async < parallel_output_complete_action > (opts.all_localities[i], opts.data_dir, this_fname, get_time(), cycle);
							}
							hpx::wait_all(futs);
							grid::output_header(opts.data_dir, silo_name, get_time(), cycle, opts.all_localities.size());
						} else {
							this_fname = silo_name + std::string(".silo");
							grid::output(
									my_list, opts.data_dir, this_fname, current_time, get_rotation_count() / opts.output_dt);
						}
					}
					printf("Done...\n");

					if( !opts.parallel_silo) {
					}
					printf("Loaded checkpoint file\n");
					my_list = decltype(my_list)();
				}

				return my_list;
			}, std::move(futs));

}
grid::output_list_type node_server::output(std::string dname, std::string fname,
		int cycle) const {
	if (is_refined) {
		std::array<hpx::future<grid::output_list_type>, NCHILD> futs;
		integer index = 0;
		for (auto i = children.begin(); i != children.end(); ++i) {
			futs[index++] = i->output(dname, fname, cycle);
		}

		auto i = futs.begin();
		grid::output_list_type my_list = i->get();
		if (opts.parallel_silo) {
			parallel_output_gather(std::move(my_list));
			for (++i; i != futs.end(); ++i) {
				i->get();
			}
		} else {
			for (++i; i != futs.end(); ++i) {
				grid::merge_output_lists(my_list, i->get());
			}
		}
		if (my_location.level() == 0) {
			std::string this_fname;
			printf("Outputing...\n");
			if (opts.parallel_silo) {
				std::string this_dname = dname + fname
						+ std::string(".silo.data/");
				//printf("node_server::output (mkdir): this_dname('%s')\n", this_dname.c_str());
				if (system(
						(std::string("mkdir -p ") + this_dname
								+ std::string("\n")).c_str()) != 0) {
					abort_error()
					;
				}
				const auto sz = opts.all_localities.size();
				std::vector<hpx::future<void>> futs(sz);
				for (integer i = 0; i != sz; ++i) {
					this_fname = fname + std::string(".") + std::to_string(i)
							+ std::string(".silo");
					futs[i] =
							hpx::async < parallel_output_complete_action
									> (opts.all_localities[i], this_dname, this_fname, get_time(), cycle);
				}
				hpx::wait_all(futs);
				grid::output_header(this_dname, fname, get_time(), cycle,
						opts.all_localities.size());
			} else {
				this_fname = fname + std::string(".silo");
				grid::output(my_list, dname, this_fname, get_time(), cycle);
			}
			printf("Done...\n");
		}
		return my_list;
	} else {
		auto l = grid_ptr->get_output_list();
		if (opts.parallel_silo) {
			parallel_output_gather(std::move(l));
		}
		return l;
	}
}

integer node_server::regrid_gather(bool rebalance_only) {
	integer count = integer(1);

	if (is_refined) {
		if (!rebalance_only) {
			/* Turning refinement off */
			if (refinement_flag == 0) {
				std::fill_n(children.begin(), NCHILD, node_client());
				is_refined = false;
			}
		}

		if (is_refined) {
			std::array<hpx::future<integer>, NCHILD> futs;
			integer index = 0;
			for (auto& child : children) {
				futs[index++] = child.regrid_gather(rebalance_only);
			}
			auto futi = futs.begin();
			for (auto const& ci : geo::octant::full_set()) {
				auto child_cnt = futi->get();
				++futi;
				child_descendant_count[ci] = child_cnt;
				count += child_cnt;
			}
		} else {
			for (auto const& ci : geo::octant::full_set()) {
				child_descendant_count[ci] = 0;
			}
		}
	} else if (!rebalance_only) {
		//		if (grid_ptr->refine_me(my_location.level())) {
		if (refinement_flag != 0) {
			refinement_flag = 0;
			count += NCHILD;

			/* Turning refinement on*/
			is_refined = true;

			for (auto& ci : geo::octant::full_set()) {
				child_descendant_count[ci] = 1;

//                 children[ci] = create_child(hpx::find_here(), ci);

			}
		}
	}

	return count;
}

hpx::future<hpx::id_type> node_server::create_child(
		hpx::id_type const& locality, integer ci) {
	return hpx::new_ < node_server
			> (hpx::find_here(), my_location.get_child(ci), me, current_time, step_num, hcycle).then(
					[this, ci](hpx::future<hpx::id_type>&& child_idf)
					{
						hpx::id_type child_id = child_idf.get();
						node_client child = child_id;
						{
							std::array<integer, NDIM> lb = {2 * BW, 2 * BW, 2 * BW};
							std::array<integer, NDIM> ub;
							lb[XDIM] += (1 & (ci >> 0)) * (INX);
							lb[YDIM] += (1 & (ci >> 1)) * (INX);
							lb[ZDIM] += (1 & (ci >> 2)) * (INX);
							for (integer d = 0; d != NDIM; ++d) {
								ub[d] = lb[d] + (INX);
							}
							if (current_time > ZERO)
							{
								auto prolong = grid_ptr->get_prolong(lb, ub);
								child.set_grid(std::move(prolong)).get();
							}
						}
						return child_id;
					});
}

hpx::future<void> node_server::regrid_scatter(integer a_, integer total) {
	refinement_flag = 0;
	std::array<hpx::future<void>, geo::octant::count()> futs;
	if (is_refined) {
		integer a = a_;
		++a;
		integer index = 0;
		for (auto& ci : geo::octant::full_set()) {
			const integer loc_index = a * options::all_localities.size()
					/ total;
			const auto child_loc = options::all_localities[loc_index];
			if (children[ci].empty()) {
				futs[index++] = create_child(child_loc, ci).then(
						[this, ci, a, total](hpx::future<hpx::id_type>&& child)
						{
							children[ci] = child.get();
							return children[ci].regrid_scatter(a, total);
						});
			} else {
				const hpx::id_type id = children[ci].get_gid();
				integer current_child_id =
						hpx::naming::get_locality_id_from_gid(id.get_gid());
				auto current_child_loc =
						options::all_localities[current_child_id];
				if (child_loc != current_child_loc) {
					futs[index++] =
							children[ci].copy_to_locality(child_loc).then(
									[this, ci, a, total](hpx::future<hpx::id_type>&& child)
									{
										children[ci] = child.get();
										return children[ci].regrid_scatter(a, total);
									});
				} else {
					futs[index++] = children[ci].regrid_scatter(a, total);
				}
			}
			a += child_descendant_count[ci];
		}
	}
	clear_family();
	if (is_refined) {
		return hpx::when_all(futs);
	} else {
		return hpx::make_ready_future();
	}
}

integer node_server::regrid(const hpx::id_type& root_gid, bool rb) {
	hpx::util::high_resolution_timer timer;
	assert(grid_ptr != nullptr);
	printf("-----------------------------------------------\n");
	if (!rb) {
		printf("checking for refinement\n");
		check_for_refinement().get();
	}
	printf("regridding\n");
	real tstart = timer.elapsed();
	integer a = regrid_gather(rb);
	real tstop = timer.elapsed();
	printf("Regridded tree in %f seconds\n", real(tstop - tstart));
	printf("rebalancing %i nodes\n", int(a));
	tstart = timer.elapsed();
	regrid_scatter(0, a).get();
	tstop = timer.elapsed();
	printf("Rebalanced tree in %f seconds\n", real(tstop - tstart));
	assert(grid_ptr != nullptr);
	std::vector<hpx::id_type> null_neighbors(geo::direction::count());
	tstart = timer.elapsed();
	printf("forming tree connections\n");
	form_tree(hpx::unmanaged(root_gid), hpx::invalid_id, null_neighbors).get();
	tstop = timer.elapsed();
	printf("Formed tree in %f seconds\n", real(tstop - tstart));
	double elapsed = timer.elapsed();
	printf(
			"regrid done in %f seconds\n---------------------------------------\n",
			elapsed);
	return a;
}

std::map<integer, std::vector<char> > node_server::save_local(integer& cnt,
		std::string const& filename, hpx::future<void>& child_fut) const {

	std::map<integer, std::vector<char> > result;
	char flag = is_refined ? '1' : '0';
	integer my_cnt = cnt;

	// Call save on children that are non local, for all
	// locals, get the pointer
	std::vector<hpx::future<void>> child_futs;
	std::vector<hpx::util::tuple<integer, integer, hpx::future<node_server*>>>local_children;
	if (is_refined) {
		child_futs.resize(NCHILD);
		local_children.reserve(NCHILD);
		integer i = cnt + 1;
		for (auto& ci : geo::octant::full_set()) {
			if (!children[ci].is_local()) {
				child_futs[ci] = children[ci].save(i, filename);
			} else {
				local_children.emplace_back(ci, i, children[ci].get_ptr());
			}
			i += child_descendant_count[ci];
		}
	}

	std::vector<char> out_buffer;
	out_buffer.reserve(4096 * 5);
	{
		typedef hpx::util::container_device<std::vector<char> > io_device_type;
		boost::iostreams::stream<io_device_type> strm(out_buffer);

		// determine record size
		write(strm, flag);
		integer value = ++cnt;
		std::array<integer, NCHILD> values;
		for (auto const& ci : geo::octant::full_set()) {
			if (ci != 0 && is_refined) {
				value += child_descendant_count[ci - 1];
			}
			values[ci] = value;
			write(strm, value);
		}
		save_me(strm);
	}
	result.emplace(my_cnt, std::move(out_buffer));

	if (is_refined) {
		for (auto& child : local_children) {
			auto ci = hpx::util::get < 0 > (child);
			auto cnt = hpx::util::get < 1 > (child);
			auto cc = hpx::util::get < 2 > (child).get();
			auto child_result = cc->save_local(cnt, filename, child_futs[ci]);
			for (auto && d : child_result) {
				result.emplace(std::move(d));
			}
		}
		child_fut = hpx::when_all(child_futs);
	}

	return result;
}

hpx::future<void> node_server::save(integer cnt,
		std::string const& filename) const {
	// Create file and save metadata if location == 0
	if (my_location.level() == 0) {
		// run output on separate thread
		hpx::threads::run_as_os_thread([&]() {
			FILE *fp = fopen(filename.c_str(), "wb");
			fclose(fp);
		}).get();
	}

	hpx::future<void> child_fut;
	std::map<integer, std::vector<char> > result = save_local(cnt, filename,
			child_fut);

	// run output on separate thread
	auto fut =
			hpx::threads::run_as_os_thread([&]() {
				// write all of the buffers to file
					integer record_size = 0;
					FILE* fp = fopen(filename.c_str(), "rb+");
					for (auto const& d : result) {
						if (record_size == 0) {
							record_size = d.second.size();
						}
						else {
							assert(record_size == d.second.size());
						}
						fseek(fp, record_size * d.first, SEEK_SET);
						fwrite(d.second.data(), sizeof(char), d.second.size(), fp);
					}

					if (my_location.level() == 0) {
						std::size_t total = 1;
						for (auto& ci: geo::octant::full_set())
						{
							total += child_descendant_count[ci];
						}
						fseek(fp, record_size * total, SEEK_SET);
						fwrite(&record_size, sizeof(integer), 1, fp);
						fwrite(outflows.data(), sizeof(real), physics::NF, fp);
						printf("Saved %li grids to checkpoint file\n", (long int) total);
					}
					fclose(fp);
				});

	return hpx::when_all(fut, child_fut);
}

void node_server::set_aunt(const hpx::id_type& aunt, const geo::face& face) {
	aunts[face] = aunt;
}

void node_server::set_grid(const std::vector<simd_vector>& data) {
	grid_ptr->set_prolong(data);
}

hpx::future<void> node_server::check_for_refinement() {
	bool rc = false;
	std::array<hpx::future<void>, NCHILD + 1> futs;
	integer index = 0;
	if (is_refined) {
		for (auto& child : children) {
			futs[index++] = child.check_for_refinement();
		}
	}
	all_hydro_bounds();
	if (!rc) {
		rc = grid_ptr->refine_me(my_location.level());
	}
	if (rc) {
		if (refinement_flag++ == 0) {
			if (!parent.empty()) {
				futs[index++] = parent.force_nodes_to_exist(
						my_location.get_neighbors());
			}
		}
	}
	if (is_refined) {
		return hpx::when_all(futs);
	} else {
		return hpx::make_ready_future();
	}
}

hpx::future<hpx::id_type> node_server::copy_to_locality(
		const hpx::id_type& id) {

	std::vector<hpx::id_type> cids;
	if (is_refined) {
		cids.resize(NCHILD);
		for (auto& ci : geo::octant::full_set()) {
			cids[ci] = children[ci].get_gid();
		}
	}
	auto rc =
			hpx::new_ < node_server
					> (id, my_location, step_num, is_refined, current_time, child_descendant_count, std::move(
							*grid_ptr), cids, std::size_t(hcycle));
	clear_family();
	return rc;
}

void node_server::force_nodes_to_exist(std::vector<node_location>&& locs) {
	std::vector<hpx::future<void>> futs;
	std::vector<node_location> parent_list;
	std::vector<std::vector<node_location>> child_lists(NCHILD);

	futs.reserve(geo::octant::count() + 2);
	parent_list.reserve(locs.size());

	integer index = 0;
	for (auto& loc : locs) {
		assert(loc != my_location);
		if (loc.is_child_of(my_location)) {
			if (refinement_flag++ == 0 && !parent.empty()) {
				futs.push_back(
						parent.force_nodes_to_exist(
								my_location.get_neighbors()));
			}
			if (is_refined) {
				for (auto& ci : geo::octant::full_set()) {
					if (loc.is_child_of(my_location.get_child(ci))) {
						if (child_lists[ci].empty()) {
							child_lists[ci].reserve(locs.size());
						}
						child_lists[ci].push_back(loc);
						break;
					}
				}
			}
		} else {
			assert(!parent.empty());
			parent_list.push_back(loc);
		}
	}
	for (auto& ci : geo::octant::full_set()) {
		if (is_refined && child_lists[ci].size()) {
			futs.push_back(
					children[ci].force_nodes_to_exist(
							std::move(child_lists[ci])));
		}
	}
	if (parent_list.size()) {
		futs.push_back(parent.force_nodes_to_exist(std::move(parent_list)));
	}

	wait_all_and_propagate_exceptions(futs);
}

hpx::future<void> node_server::form_tree(hpx::id_type self_gid,
		hpx::id_type parent_gid, std::vector<hpx::id_type> neighbor_gids) {
	for (auto& dir : geo::direction::full_set()) {
		neighbors[dir] = std::move(neighbor_gids[dir]);
	}
	me = std::move(self_gid);
	parent = std::move(parent_gid);
	if (is_refined) {
		std::array<hpx::future<void>, NCHILD> cfuts;
		integer index = 0;
		amr_flags.resize(NCHILD);
		for (integer cx = 0; cx != 2; ++cx) {
			for (integer cy = 0; cy != 2; ++cy) {
				for (integer cz = 0; cz != 2; ++cz) {
					std::array<hpx::future<hpx::id_type>,
							geo::direction::count()> child_neighbors_f;
					const integer ci = cx + 2 * cy + 4 * cz;
					for (integer dx = -1; dx != 2; ++dx) {
						for (integer dy = -1; dy != 2; ++dy) {
							for (integer dz = -1; dz != 2; ++dz) {
								if (!(dx == 0 && dy == 0 && dz == 0)) {
									const integer x = cx + dx + 2;
									const integer y = cy + dy + 2;
									const integer z = cz + dz + 2;
									geo::direction i;
									i.set(dx, dy, dz);
									auto& ref = child_neighbors_f[i];
									auto other_child = (x % 2) + 2 * (y % 2)
											+ 4 * (z % 2);
									if (x / 2 == 1 && y / 2 == 1
											&& z / 2 == 1) {
										ref =
												hpx::make_ready_future
														< hpx::id_type
														> (hpx::unmanaged(
																children[other_child].get_gid()));
									} else {
										geo::direction dir =
												geo::direction(
														(x / 2)
																+ NDIM
																		* ((y
																				/ 2)
																				+ NDIM
																						* (z
																								/ 2)));
										node_location parent_loc =
												my_location.get_neighbor(dir);
										ref = neighbors[dir].get_child_client(
												parent_loc, other_child);
									}
								}
							}
						}
					}
					cfuts[index++] =
							hpx::dataflow(hpx::launch::sync,
									[this, ci](std::array<hpx::future<hpx::id_type>, geo::direction::count()>&& cns) {
										std::vector<hpx::id_type> child_neighbors(geo::direction::count());
										for (auto& dir : geo::direction::full_set()) {
											child_neighbors[dir] = cns[dir].get();
											amr_flags[ci][dir] = bool(child_neighbors[dir] == hpx::invalid_id);
										}
										return children[ci].form_tree(hpx::unmanaged(children[ci].get_gid()),
												me.get_gid(), std::move(child_neighbors));
									}, std::move(child_neighbors_f));
				}
			}
		}
		return hpx::when_all(cfuts);
	} else {
		std::vector<hpx::future<void>> nfuts;
		nfuts.reserve(NFACE);
		for (auto& f : geo::face::full_set()) {
			const auto& neighbor = neighbors[f.to_direction()];
			if (!neighbor.empty()) {
				nfuts.push_back(
						neighbor.set_child_aunt(me.get_gid(), f ^ 1).then(
								[this, f](hpx::future<bool>&& n)
								{
									nieces[f] = n.get();
								}));
			} else {
				nieces[f] = false;
			}
		}
		return hpx::when_all(nfuts);
	}
}

hpx::id_type node_server::get_child_client(const geo::octant& ci) {
	if (is_refined) {
		return children[ci].get_gid();
	} else {
		return hpx::invalid_id;
	}
}

bool node_server::set_child_aunt(const hpx::id_type& aunt,
		const geo::face& face) const {
	if (is_refined) {
		std::array<hpx::future<void>, geo::octant::count() / 2> futs;
		integer index = 0;
		for (auto const& ci : geo::octant::face_subset(face)) {
			futs[index++] = children[ci].set_aunt(aunt, face);
		}
		wait_all_and_propagate_exceptions(futs);
	}
	return is_refined;
}

std::uintptr_t node_server::get_ptr() {
	return reinterpret_cast<std::uintptr_t>(this);
}
