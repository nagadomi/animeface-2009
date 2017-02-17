require "mkmf"
dir_config("/usr/local/include", "/usr/local/lib")
unless ( have_header('nv_core.h') && have_library('nvxs', 'nv_matrix_alloc'))
  $stderr.puts("error: Can't locate libnvxs")
  exit(-1)
end

create_makefile("AnimeFace")
