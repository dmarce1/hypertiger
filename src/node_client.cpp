/*
 * node_client.cpp
 *
 *  Created on: Jun 12, 2015
 *      Author: dmarce1
 */

#include "node_server.hpp"
#include "options.hpp"


#include <hpx/lcos/broadcast.hpp>
#include <hpx/runtime/get_colocation_id.hpp>

hpx::id_type node_client::get_gid() const {
	return id;
}

hpx::id_type node_client::get_unmanaged_gid() const {
	return unmanaged;
}

node_client& node_client::operator=(hpx::future<hpx::id_type>&& fut) {
	id = fut.get();
	if (!empty()) {
		unmanaged = hpx::id_type(id.get_gid(), hpx::id_type::unmanaged);
	}
	return *this;
}

node_client& node_client::operator=(const hpx::id_type& _id) {
	id = _id;
	if (!empty()) {
		unmanaged = hpx::id_type(id.get_gid(), hpx::id_type::unmanaged);
	}
	return *this;
}

bool node_client::is_local() const {
	return hpx::get_colocation_id(id).get() == hpx::find_here();
}

node_client::node_client(hpx::future<hpx::id_type>&& fut) {
	id = fut.get();
	if (!empty()) {
		unmanaged = hpx::id_type(id.get_gid(), hpx::id_type::unmanaged);
	}
}

node_client::node_client(const hpx::id_type& _id) {
	id = _id;
	if (!empty()) {
		unmanaged = hpx::id_type(id.get_gid(), hpx::id_type::unmanaged);
	}
}

node_client::node_client() {
}

bool node_client::empty() const {
	return get_gid() == hpx::invalid_id;
}

hpx::future<grid::output_list_type> node_client::load(integer i, integer total,
		integer rec_size, bool do_o, std::string s) const {
	return hpx::async<typename node_server::load_action>(get_unmanaged_gid(), i,
			total, rec_size, do_o, s);
}

hpx::future<grid::output_list_type> node_client::output(std::string dname,
		std::string fname, int cycle) const {
	return hpx::async<typename node_server::output_action>(get_unmanaged_gid(),
			dname, fname, cycle);
}

hpx::future<integer> node_client::regrid_gather(bool rb) const {
	return hpx::async<typename node_server::regrid_gather_action>(
			get_unmanaged_gid(), rb);
}

hpx::future<void> node_client::regrid_scatter(integer a, integer b) const {
	return hpx::async<typename node_server::regrid_scatter_action>(
			get_unmanaged_gid(), a, b);
}

hpx::future<void> node_client::save(integer i, std::string s) const {
	return hpx::async<typename node_server::save_action>(get_unmanaged_gid(), i,
			s);
}

hpx::future<void> node_client::set_aunt(const hpx::id_type& aunt,
		const geo::face& f) const {
	return hpx::async<typename node_server::set_aunt_action>(
			get_unmanaged_gid(), aunt, f);
}

hpx::future<void> node_client::set_grid(std::vector<simd_vector>&& g) const {
	return hpx::async<typename node_server::set_grid_action>(
			get_unmanaged_gid(), std::move(g));
}

hpx::future<void> node_client::check_for_refinement() const {
	return hpx::async<typename node_server::check_for_refinement_action>(
			get_unmanaged_gid());
}

hpx::future<hpx::id_type> node_client::copy_to_locality(
		const hpx::id_type& id) const {
	return hpx::async<typename node_server::copy_to_locality_action>(get_gid(),
			id);
}

hpx::future<void> node_client::force_nodes_to_exist(
		std::vector<node_location>&& locs) const {
	return hpx::async<typename node_server::force_nodes_to_exist_action>(
			get_unmanaged_gid(), std::move(locs));
}

hpx::future<void> node_client::form_tree(hpx::id_type&& id1, hpx::id_type&& id2,
		std::vector<hpx::id_type>&& ids) {
	return hpx::async<typename node_server::form_tree_action>(
			get_unmanaged_gid(), std::move(id1), std::move(id2), std::move(ids));
}

hpx::future<hpx::id_type> node_client::get_child_client(
		const node_location& parent_loc, const geo::octant& ci) {
	hpx::future < hpx::id_type > rfut;
	if (get_gid() != hpx::invalid_id) {
		rfut = hpx::async<typename node_server::get_child_client_action>(
				get_unmanaged_gid(), ci);
		;
	} else {
		auto tmp = hpx::invalid_id;
		rfut = hpx::make_ready_future < hpx::id_type > (std::move(tmp));
	}
	return rfut;
}

hpx::future<bool> node_client::set_child_aunt(const hpx::id_type& aunt,
		const geo::face& f) const {
	return hpx::async<typename node_server::set_child_aunt_action>(
			get_unmanaged_gid(), aunt, f);
}

hpx::future<node_server*> node_client::get_ptr() const {
	return hpx::async<typename node_server::get_ptr_action>(get_unmanaged_gid()).then(
			[](hpx::future<std::uintptr_t>&& fut) {
				return reinterpret_cast<node_server*>(fut.get());
			});
}

void node_client::send_hydro_boundary(std::vector<simd_vector>&& data,
		const geo::direction& dir, std::size_t cycle) const {
	hpx::apply<typename node_server::send_hydro_boundary_action>(
			get_unmanaged_gid(), std::move(data), dir, cycle);
}

void node_client::send_hydro_children(std::vector<real>&& data,
		const geo::octant& ci, std::size_t cycle) const {
	hpx::apply<typename node_server::send_hydro_children_action>(
			get_unmanaged_gid(), std::move(data), ci, cycle);
}

void node_client::send_hydro_flux_correct(std::vector<real>&& data,
		const geo::face& face, const geo::octant& ci) const {
	hpx::apply<typename node_server::send_hydro_flux_correct_action>(
			get_unmanaged_gid(), std::move(data), face, ci);
}

hpx::future<std::pair<std::vector<simd_vector>, real>> node_client::step(
		integer steps) const {
	return hpx::async<typename node_server::step_action>(get_unmanaged_gid(),
			steps);
}

void node_client::timestep_driver_ascend(real dt) const {
	hpx::apply<typename node_server::timestep_driver_ascend_action>(
			get_unmanaged_gid(), dt);
}

void node_client::set_local_timestep(integer idx, real dt) const {
	hpx::apply<typename node_server::set_local_timestep_action>(
			get_unmanaged_gid(), idx, dt);
}

