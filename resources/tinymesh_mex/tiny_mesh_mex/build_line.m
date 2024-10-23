
mex -v taubin_smooth_tiny_mesh_mex.cpp -IF:\EM_astrocyte\codes\EM_dendrite\resources\tinymesh\src\tinymesh -L"F:\EM_astrocyte\codes\EM_dendrite\resources\tinymesh\build\lib" -lF:\EM_astrocyte\codes\EM_dendrite\resources\tinymesh\build\lib\tinymesh.lib 



mex -v taubin_smooth_tiny_mesh_mex.cpp -I'F:\EM_astrocyte\codes\EM_dendrite\resources\tinymesh\src\tinymesh'


mex CXXFLAGS='$CXXFLAGS -std=c++14' taubin_smooth_tiny_mesh_mex.cpp -I/work/boyu/EM_astrocyte/EM_dendrite/resources/tiny_mesh_ubuntu/tinymesh/tinymesh/src/tinymesh ...
    -I/work/boyu/EM_astrocyte/EM_dendrite/resources/tiny_mesh_ubuntu/tinymesh/tinymesh/src/tinymesh/ext/eigen ...
    -L/work/boyu/EM_astrocyte/EM_dendrite/resources/tiny_mesh_ubuntu/tinymesh/tinymesh/build/lib ...
    -ltinymesh


