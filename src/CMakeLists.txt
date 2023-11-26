
set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(APPEND srcs
   "${dir}/lwoniom_parse.F90"
   "${dir}/lwoniom_structure.f90"
   "${dir}/lwoniom_covrad.f90"
   "${dir}/lwoniom_engrad.f90"
   "${dir}/lwoniom_hessian.f90"
   "${dir}/lwoniom_setup.f90"
   "${dir}/lwoniom_interface.f90"
)

set(srcs ${srcs} PARENT_SCOPE)
