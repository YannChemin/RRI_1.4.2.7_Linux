# RRI_1.4.2.7_Linux

To compile as is:

~~~shell
make clean
make mod
make
~~~

To compile as debug (full terminal info):

~~~shell
cp debug/RRI.f90 .
cp debug/RRI_Read.f90 .
make clean
make mod
make
~~~

To compile as operational (no terminal info):

~~~shell
cp operational/RRI.f90 .
cp operational/RRI_Read.f90 .
make clean
make mod
make
~~~

To run

~~~shell
mv 0_rri_1_4_2_7_Linux solo30s/
cd solo30s
chmod 755 0_rri_1_4_2_7_Linux
./0_rri_1_4_2_7_Linux
~~~

or compile + run

~~~shell
bash ./compile.sh
~~~

The same type of modifications can be done to the type of parallelisation

To compile with standard OpenMP:

~~~shell
cp openmp/RRI_GW.f90 .
cp openmp/RRI_Riv.f90 .
cp openmp/RRI_Slope.f90 .
make clean
make mod
make
~~~

To compile with OpenMP offloading to GPU (experimental):

~~~shell
cp gpuoffloading/RRI_GW.f90 .
cp gpuoffloading/RRI_Riv.f90 .
cp gpuoffloading/RRI_Slope.f90 .
make clean
make mod
make
~~~
