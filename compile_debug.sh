# Replace code with DEBUG version
cp debug/RRI_Read.f90 RRI_Read.f90
cp debug/RRI.f90 RRI.f90
# Compile
make clean
make mod
make
# Move Exec
mv 0_rri_1_4_2_7_Linux solo30s/
# Go to location
cd solo30s
# Make Exec permisions
chmod 755 0_rri_1_4_2_7_Linux
# Run the model
time ./0_rri_1_4_2_7_Linux
