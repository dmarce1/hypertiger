/*
 * node_client.hpp
 *
 *  Created on: Jun 11, 2015
 *      Author: dmarce1
 */

#ifndef NODE_CLIENT_HPP_
#define NODE_CLIENT_HPP_

#include "defs.hpp"
#include "node_location.hpp"
#include "grid.hpp"
#include "geometry.hpp"

class node_server;

namespace hpx {
using mutex = hpx::lcos::local::spinlock;
}

class node_client {

private:
	hpx::id_type id;
	hpx::id_type unmanaged;
public:
	bool is_local() const;
	template<class Arc>
	void load(Arc& arc, unsigned) {
		arc & id;
		if (!empty()) {
			unmanaged = hpx::id_type(id.get_gid(), hpx::id_type::unmanaged);
		}
	}

	template<class Arc>
	void save(Arc& arc, unsigned) const {
		arc & id;
	}
	HPX_SERIALIZATION_SPLIT_MEMBER();

	bool empty() const;
	hpx::id_type get_gid() const;
	hpx::id_type get_unmanaged_gid() const;
	node_client& operator=(hpx::future<hpx::id_type>&& fut);
	node_client& operator=(const hpx::id_type& _id);
	node_client(hpx::future<hpx::id_type>&& fut);
	node_client(const hpx::id_type& _id);
	void send_hydro_children(std::vector<real>&&, const geo::octant& ci,
			std::size_t cycle) const;
	void send_hydro_flux_correct(std::vector<real>&&, const geo::face& face,
			const geo::octant& ci) const;
	hpx::future<grid::output_list_type> load(integer, integer, integer,
			bool do_output, std::string) const;
	hpx::future<grid::output_list_type> output(std::string dname,
			std::string fname, int) const;
	node_client();
	hpx::future<bool> set_child_aunt(const hpx::id_type&,
			const geo::face&) const;
	hpx::future<void> set_aunt(const hpx::id_type&, const geo::face&) const;
	hpx::future<node_server*> get_ptr() const;
	hpx::future<void> form_tree(hpx::id_type&&, hpx::id_type&&,
			std::vector<hpx::id_type>&&);
	hpx::future<hpx::id_type> get_child_client(const node_location& parent_loc,
			const geo::octant&);
	hpx::future<void> regrid_scatter(integer, integer) const;
	hpx::future<integer> regrid_gather(bool) const;
	void send_hydro_boundary(std::vector<simd_vector>&&, const geo::direction& dir,
			std::size_t cycle) const;
	hpx::future<std::pair<std::vector<simd_vector>, real>> step(integer) const;
	hpx::future<hpx::id_type> copy_to_locality(const hpx::id_type&) const;
	hpx::future<void> set_grid(std::vector<simd_vector>&&) const;
	void timestep_driver_ascend(real) const;
	void set_local_timestep(integer, real) const;
	hpx::future<grid::output_list_type> output() const;
	hpx::future<void> save(integer, std::string) const;
	hpx::future<void> check_for_refinement() const;
	hpx::future<void> force_nodes_to_exist(
			std::vector<node_location>&& loc) const;

	void report_timing() const;

};
#endif /* NODE_CLIENT_HPP_ */
