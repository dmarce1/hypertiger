/*
 * node_server_actions.cpp
 *
 *  Created on: Apr 15, 2017
 *      Author: dminacore
 */


#include "node_server.hpp"

HPX_REGISTER_COMPONENT(hpx::components::managed_component<node_server>, node_server);

typedef node_server::set_child_aunt_action set_child_aunt_action_type;
typedef node_server::get_ptr_action get_ptr_action_type;
typedef node_server::send_hydro_boundary_action send_hydro_boundary_action_type;
typedef node_server::send_hydro_children_action send_hydro_children_action_type;
typedef node_server::send_hydro_flux_correct_action send_hydro_flux_correct_action_type;
typedef node_server::step_action step_action_type;
typedef node_server::timestep_driver_ascend_action timestep_driver_ascend_action_type;
typedef node_server::set_local_timestep_action set_local_timestep_action_type;
typedef node_server::load_action load_action_type;
typedef node_server::output_action output_action_type;
typedef node_server::regrid_gather_action regrid_gather_action_type;
typedef node_server::regrid_scatter_action regrid_scatter_action_type;
typedef node_server::save_action save_action_type;
typedef node_server::set_aunt_action set_aunt_action_type;
typedef node_server::set_grid_action set_grid_action_type;
typedef node_server::check_for_refinement_action check_for_refinement_action_type;
typedef node_server::copy_to_locality_action copy_to_locality_action_type;
typedef node_server::force_nodes_to_exist_action force_nodes_to_exist_action_type;
typedef node_server::form_tree_action form_tree_action_type;
typedef node_server::get_child_client_action get_child_client_action_type;

HPX_REGISTER_ACTION (set_child_aunt_action_type);
HPX_REGISTER_ACTION (get_ptr_action_type);
HPX_REGISTER_ACTION (send_hydro_boundary_action_type);
HPX_REGISTER_ACTION (send_hydro_children_action_type);
HPX_REGISTER_ACTION (send_hydro_flux_correct_action_type);
HPX_REGISTER_ACTION (step_action_type);
HPX_REGISTER_ACTION (timestep_driver_ascend_action_type);
HPX_REGISTER_ACTION (set_local_timestep_action_type);
HPX_REGISTER_ACTION (load_action_type);
HPX_REGISTER_ACTION (output_action_type);
HPX_REGISTER_ACTION (regrid_gather_action_type);
HPX_REGISTER_ACTION (regrid_scatter_action_type);
HPX_REGISTER_ACTION (save_action_type);
HPX_REGISTER_ACTION (set_aunt_action_type);
HPX_REGISTER_ACTION (set_grid_action_type);
HPX_REGISTER_ACTION (check_for_refinement_action_type);
HPX_REGISTER_ACTION (copy_to_locality_action_type);
HPX_REGISTER_ACTION (force_nodes_to_exist_action_type);
HPX_REGISTER_ACTION (form_tree_action_type);
HPX_REGISTER_ACTION (get_child_client_action_type);
