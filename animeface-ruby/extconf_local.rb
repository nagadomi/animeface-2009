require "mkmf"

dir_config("nvxs", "../install/include", "../install/lib")
rpath=File.expand_path(File.join(File.dirname(__FILE__), '../install/lib'))
$LDFLAGS << " -Wl,-rpath,#{rpath}"
unless ( have_header('nv_core.h') && have_library('nvxs', 'nv_matrix_alloc'))
  $stderr.puts("error: Can't locate libnvxs")
  exit(-1)
end
create_makefile("AnimeFace")
