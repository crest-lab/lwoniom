
set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(APPEND srcs3
  "${dir}/example.f90"
  "${dir}/subprocess.f90"
  "${dir}/main.f90"
)

set(srcs3 ${srcs3} PARENT_SCOPE)
