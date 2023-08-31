# This file is part of lwoniom.
# SPDX-Identifier: LGPL-3.0-or-later
#
# lwoniom is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# lwoniom is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with lwoniom.  If not, see <https://www.gnu.org/licenses/>.

set(_lib "toml-f")
set(_pkg "TOML-F")
set(_url "https://github.com/toml-f/toml-f")

if(NOT DEFINED "${_pkg}_FIND_METHOD")
  if(DEFINED "${PROJECT_NAME}-dependency-method")
    set("${_pkg}_FIND_METHOD" "${${PROJECT_NAME}-dependency-method}")
  else()
    set("${_pkg}_FIND_METHOD" "cmake" "pkgconf" "subproject" "fetch")
  endif()
  set("_${_pkg}_FIND_METHOD")
endif()

include("${CMAKE_CURRENT_LIST_DIR}/lwoniom-utils.cmake")

lwoniom_find_package("${_lib}" "${${_pkg}_FIND_METHOD}" "${_url}")

if(TARGET "toml-f::toml-f")
  set (found TRUE)
else()
  set (found FALSE)
endif()
message("-- Found toml-f: ${found}")

if(DEFINED "_${_pkg}_FIND_METHOD")
  unset("${_pkg}_FIND_METHOD")
  unset("_${_pkg}_FIND_METHOD")
endif()
unset(_lib)
unset(_pkg)
unset(_url)
