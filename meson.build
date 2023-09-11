project(
  'lwoniom',
  'fortran', 'c',
  version: '0.0.1',
  license: 'LGPL-3.0-or-later',
  meson_version: '>=0.63',
  default_options: [
    'buildtype=debugoptimized',
    'default_library=static',
    'c_link_args=-static',
    'fortran_link_args=-static',
  ],
)
install = not (meson.is_subproject() and get_option('default_library') == 'static')

# General configuration information
exe_deps = []
subdir('config')

# Documentation
#subdir('docs')

# Collect source of the project
prog = []
tester = []
srcs = []
subdir('src')


# Create library target
lwoniom_lib = library(
  meson.project_name(),
  sources: srcs,
  dependencies: exe_deps,
  include_directories: include_directories('include'),
)

# Export as dependency
lwoniom_inc = [include_directories('include'), lwoniom_lib.private_dir_include()]
lwoniom_dep = declare_dependency(
  link_with: lwoniom_lib,
  include_directories: lwoniom_inc,
  dependencies: exe_deps,
  variables: {'includedir': meson.current_source_dir() / 'include'},
)


# Create tester target
if get_option('build_test')
subdir('testprog')
lwoniom_exe = executable(
  'lwoniom-tester',
  sources: tester,
  dependencies: lwoniom_dep,
  install: install,
  link_language: 'fortran',
)
endif

if get_option('build_exe')
subdir('app')
lwoniom_exe = executable(
  'lwoniom-app',
  sources: prog,
  dependencies: lwoniom_dep,
  install: install,
  link_language: 'fortran',
)
endif


