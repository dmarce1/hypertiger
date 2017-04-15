/*
 * node_server.hpp
 *
 *  Created on: Jun 11, 2015
 *      Author: dmarce1
 */

#ifndef NODE_SERVER_HPP_
#define NODE_SERVER_HPP_

#include "defs.hpp"
#include "node_location.hpp"
#include "node_client.hpp"
#include "grid.hpp"
#include "geometry.hpp"
#include "channel.hpp"
#include "future.hpp"

#include <array>
#include <atomic>
#include <iostream>
#include <map>
#include <vector>

#include <hpx/include/components.hpp>
#include <hpx/include/serialization.hpp>

class node_server: public hpx::components::managed_component_base<node_server> {
private:
	struct sibling_hydro_type {
		std::vector<simd_vector> data;
		geo::direction direction;
	};
	std::atomic<integer> refinement_flag;
	static std::vector<real> outflows;
	static std::stack<grid::output_list_type> pending_output;
	node_location my_location;
	integer step_num;
	std::size_t hcycle;
	real current_time;
	std::shared_ptr<grid> grid_ptr; //
	bool is_refined;
	std::array<integer, NVERTEX> child_descendant_count;
	std::array<real, NDIM> xmin;
	real dx;

	node_client me;
	node_client parent;
	std::vector<node_client> neighbors;
	std::array<node_client, NCHILD> children;
	std::vector<bool> nieces;
	std::vector<node_client> aunts;

	std::vector<std::array<bool, geo::direction::count()>> amr_flags;
	hpx::lcos::local::spinlock mtx;
	std::array<channel<std::vector<real>>, NCHILD> child_hydro_channels;
	std::array<channel<sibling_hydro_type>, geo::direction::count()> sibling_hydro_channels;
	std::array<std::array<channel<std::vector<real>>, 4>, NFACE> niece_hydro_channels;
	channel<real> global_timestep_channel;
	std::array<channel<real>, NCHILD + 1> local_timestep_channels;
	hpx::mutex load_mutex;
	real dt_;

public:
	real get_time() const {
		return current_time;
	}
	real get_rotation_count() const;
	std::size_t load_me(FILE *fp, bool old_format);
	std::size_t save_me(std::ostream& strm) const;
private:

	void initialize(real);
	void send_hydro_amr_boundaries();
	void collect_hydro_boundaries();
	void exchange_interlevel_hydro_data();
	static void static_initialize();
	void all_hydro_bounds();
	void clear_family();
	hpx::future<void> exchange_flux_corrections();

	hpx::future<std::vector<simd_vector>> grid_step();

	hpx::future<std::pair<std::vector<simd_vector>, real>> local_step(integer steps);

	std::map<integer, std::vector<char> > save_local(integer& cnt,
			std::string const& filename, hpx::future<void>& child_fut) const;

public:
	static bool child_is_on_face(integer ci, integer face) {
		return (((ci >> (face / 2)) & 1) == (face & 1));
	}

	std::vector<hpx::future<void>> set_nieces_amr(const geo::face&) const;
	node_server() {
		initialize(ZERO);
	}
	~node_server() {
	}
	node_server(const node_server& other);
	node_server(const node_location&, const node_client& parent_id, real,
			std::size_t, std::size_t);
	node_server(const node_location&, integer, bool, real,
			const std::array<integer, NCHILD>&, grid,
			const std::vector<hpx::id_type>&, std::size_t);
	node_server(node_server&& other) = default;

	void save_to_file(const std::string&, std::string const& data_dir) const;
	void load_from_file(const std::string&, std::string const& data_dir);
	void load_from_file_and_output(const std::string&, const std::string&,
			std::string const& data_dir);

	grid::output_list_type output(std::string dname, std::string fname,
			int cycle) const; //
	HPX_DEFINE_COMPONENT_ACTION(node_server, output, output_action);

	static void parallel_output_gather(grid::output_list_type&&);
	static void parallel_output_complete(std::string dname, std::string fname,
			real tm, int cycle);

	integer regrid_gather(bool rebalance_only); //
	HPX_DEFINE_COMPONENT_ACTION(node_server, regrid_gather, regrid_gather_action);

	hpx::future<hpx::id_type> create_child(hpx::id_type const& locality,
			integer ci);

	hpx::future<void> regrid_scatter(integer, integer); //
	HPX_DEFINE_COMPONENT_ACTION(node_server, regrid_scatter, regrid_scatter_action);

	void recv_hydro_boundary(std::vector<simd_vector>&&, const geo::direction&,
			std::size_t cycle); //
	HPX_DEFINE_COMPONENT_DIRECT_ACTION(node_server, recv_hydro_boundary, send_hydro_boundary_action);

	void recv_hydro_children(std::vector<real>&&, const geo::octant& ci,
			std::size_t cycle); //
	HPX_DEFINE_COMPONENT_DIRECT_ACTION(node_server, recv_hydro_children, send_hydro_children_action);

	void recv_hydro_flux_correct(std::vector<real>&&, const geo::face& face,
			const geo::octant& ci); //
	HPX_DEFINE_COMPONENT_DIRECT_ACTION(node_server, recv_hydro_flux_correct, send_hydro_flux_correct_action);

	hpx::future<std::pair<std::vector<simd_vector>, real>> step(integer steps); //
	HPX_DEFINE_COMPONENT_ACTION(node_server, step, step_action);

	integer regrid(const hpx::id_type& root_gid, bool rb);

	void amr_driver();

	void set_grid(const std::vector<simd_vector>&); //
	HPX_DEFINE_COMPONENT_DIRECT_ACTION(node_server, set_grid, set_grid_action);

	hpx::future<real> timestep_driver_descend();

	void set_local_timestep(integer i, real dt); //
	HPX_DEFINE_COMPONENT_DIRECT_ACTION(node_server, set_local_timestep, set_local_timestep_action);

	void timestep_driver_ascend(real); //
	HPX_DEFINE_COMPONENT_DIRECT_ACTION(node_server, timestep_driver_ascend, timestep_driver_ascend_action);

	hpx::future<hpx::id_type> copy_to_locality(const hpx::id_type&); //
	HPX_DEFINE_COMPONENT_ACTION(node_server, copy_to_locality, copy_to_locality_action);

	hpx::id_type get_child_client(const geo::octant&); //
	HPX_DEFINE_COMPONENT_DIRECT_ACTION(node_server, get_child_client, get_child_client_action);

	hpx::future<void> form_tree(hpx::id_type, hpx::id_type,
			std::vector<hpx::id_type>); //
	HPX_DEFINE_COMPONENT_ACTION(node_server, form_tree, form_tree_action);

	std::uintptr_t get_ptr(); //
	HPX_DEFINE_COMPONENT_DIRECT_ACTION(node_server, get_ptr, get_ptr_action);

	hpx::future<grid::output_list_type> load(integer, integer, integer,
			bool do_output, std::string); //
	HPX_DEFINE_COMPONENT_ACTION(node_server, load, load_action);

	hpx::future<void> save(integer, std::string const&) const; //
	HPX_DEFINE_COMPONENT_ACTION(node_server, save, save_action);

	void set_aunt(const hpx::id_type&, const geo::face& face); //
	HPX_DEFINE_COMPONENT_DIRECT_ACTION(node_server, set_aunt, set_aunt_action);

	bool set_child_aunt(const hpx::id_type&, const geo::face& face) const; //
	HPX_DEFINE_COMPONENT_DIRECT_ACTION(node_server, set_child_aunt, set_child_aunt_action);

	hpx::future<void> check_for_refinement(); //
	HPX_DEFINE_COMPONENT_ACTION(node_server, check_for_refinement, check_for_refinement_action);

	void force_nodes_to_exist(std::vector<node_location>&& loc);HPX_DEFINE_COMPONENT_ACTION(node_server, force_nodes_to_exist, force_nodes_to_exist_action);

};

HPX_REGISTER_ACTION_DECLARATION (node_server::output_action);
HPX_REGISTER_ACTION_DECLARATION (node_server::set_grid_action);
HPX_REGISTER_ACTION_DECLARATION (node_server::force_nodes_to_exist_action);
HPX_REGISTER_ACTION_DECLARATION (node_server::check_for_refinement_action);
HPX_REGISTER_ACTION_DECLARATION (node_server::set_aunt_action);
HPX_REGISTER_ACTION_DECLARATION (node_server::set_child_aunt_action);
HPX_REGISTER_ACTION_DECLARATION (node_server::load_action);
HPX_REGISTER_ACTION_DECLARATION (node_server::save_action);
HPX_REGISTER_ACTION_DECLARATION (node_server::send_hydro_children_action);
HPX_REGISTER_ACTION_DECLARATION (node_server::send_hydro_flux_correct_action);
HPX_REGISTER_ACTION_DECLARATION (node_server::regrid_gather_action);
HPX_REGISTER_ACTION_DECLARATION (node_server::regrid_scatter_action);
HPX_REGISTER_ACTION_DECLARATION (node_server::send_hydro_boundary_action);
HPX_REGISTER_ACTION_DECLARATION (node_server::step_action);
HPX_REGISTER_ACTION_DECLARATION (node_server::copy_to_locality_action);
HPX_REGISTER_ACTION_DECLARATION (node_server::get_child_client_action);
HPX_REGISTER_ACTION_DECLARATION (node_server::form_tree_action);
HPX_REGISTER_ACTION_DECLARATION (node_server::get_ptr_action);
HPX_REGISTER_ACTION_DECLARATION (node_server::timestep_driver_ascend_action);

#endif /* NODE_SERVER_HPP_ */
