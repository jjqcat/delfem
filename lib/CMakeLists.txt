set(src_com "../src/com")
set(src_cad "../src/cad")
set(src_msh "../src/msh")
set(src_femfield "../src/femfield")
set(src_matvec "../src/matvec")
set(src_ls "../src/ls")
set(src_femls "../src/femls")
set(src_femeqn "../src/femeqn")
add_library(delfemlib STATIC 
${src_com}/vector2d.cpp 
${src_com}/vector3d.cpp 
${src_com}/drawer.cpp 
${src_com}/drawer_gl_utility.cpp 
${src_com}/octtree3d.cpp 
${src_com}/quaternion.cpp 
${src_com}/spatial_hash_grid2d.cpp 
${src_com}/spatial_hash_grid3d.cpp 
${src_com}/tri_ary_topology.cpp 
${src_com}/uglyfont.cpp 

${src_cad}/brep.cpp 
${src_cad}/brep2d.cpp 
${src_cad}/cad_edge2d_polyline.cpp 
${src_cad}/cad_elem2d.cpp 
${src_cad}/cad_elem3d.cpp 
${src_cad}/cad_obj2d.cpp 
${src_cad}/cad_obj2d_move.cpp 
${src_cad}/cad_obj3d.cpp 
${src_cad}/cad_svg.cpp 
${src_cad}/drawer_cad.cpp 
${src_cad}/drawer_cad3d.cpp

${src_msh}/drawer_msh.cpp
${src_msh}/mesh3d.cpp
${src_msh}/mesh3d_extrude.cpp
${src_msh}/mesher2d.cpp
${src_msh}/mesher2d_edit.cpp
${src_msh}/mesher3d.cpp
${src_msh}/meshkernel2d.cpp
${src_msh}/meshkernel3d.cpp
${src_msh}/surface_mesh_reader.cpp

${src_femfield}/drawer_field.cpp
${src_femfield}/drawer_field_edge.cpp
${src_femfield}/drawer_field_face.cpp
${src_femfield}/drawer_field_image_based_flow_vis.cpp
${src_femfield}/drawer_field_streamline.cpp
${src_femfield}/drawer_field_vector.cpp
${src_femfield}/elem_ary.cpp
${src_femfield}/eval.cpp
${src_femfield}/field.cpp
${src_femfield}/field_value_setter.cpp
${src_femfield}/field_world.cpp
${src_femfield}/node_ary.cpp 

${src_matvec}/mat_blkcrs.cpp 
${src_matvec}/matdia_blkcrs.cpp 
${src_matvec}/matdiafrac_blkcrs.cpp
${src_matvec}/matdiainv_blkdia.cpp
${src_matvec}/matfrac_blkcrs.cpp
${src_matvec}/matprolong_blkcrs.cpp
${src_matvec}/ordering_blk.cpp
${src_matvec}/solver_mat_iter.cpp
${src_matvec}/solver_mg.cpp
${src_matvec}/vector_blk.cpp
${src_matvec}/zmat_blkcrs.cpp
${src_matvec}/zmatdia_blkcrs.cpp
${src_matvec}/zmatdiafrac_blkcrs.cpp
${src_matvec}/zsolver_mat_iter.cpp
${src_matvec}/zvector_blk.cpp

${src_ls}/eigen_lanczos.cpp
${src_ls}/linearsystem.cpp
${src_ls}/preconditioner.cpp
${src_ls}/solver_ls_iter.cpp

${src_femls}/linearsystem_field.cpp
${src_femls}/linearsystem_fieldsave.cpp
${src_femls}/zlinearsystem.cpp
${src_femls}/zsolver_ls_iter.cpp

${src_femeqn}/eqn_advection_diffusion.cpp
${src_femeqn}/eqn_diffusion.cpp
${src_femeqn}/eqn_helmholtz.cpp
${src_femeqn}/eqn_hyper.cpp
${src_femeqn}/eqn_linear_solid2d.cpp
${src_femeqn}/eqn_linear_solid3d.cpp
${src_femeqn}/eqn_navier_stokes.cpp
${src_femeqn}/eqn_poisson.cpp
${src_femeqn}/eqn_st_venant.cpp
${src_femeqn}/eqn_stokes.cpp
${src_femeqn}/eqnsys.cpp
${src_femeqn}/eqnsys_fluid.cpp
${src_femeqn}/eqnsys_scalar.cpp
${src_femeqn}/eqnsys_shell.cpp
${src_femeqn}/eqnsys_solid.cpp
${src_femeqn}/ker_emat_tri.cpp)

