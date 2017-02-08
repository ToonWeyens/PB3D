using HDF5

print("Usage: julia concatenate_h5.jl [file_list] [output_file.h5] [dim]\n")
print("       with [file_list]      a plain file containing all the hdf5 files to merge\n")
print("            [output_file.h5] the name of the output file\n")
print("            [dim]            the dimension in which to concatenate\n")
print("\n")
print("This does NOT work when [output_file.h5] already exists\n")
print("\n")
print("Example: julia concatenate_h5.jl filelist.txt 1_sol_pert.h5 2\n")
print("         (to concatenate in the toroidal direction)\n")
print("\n")

inputfilepath=ARGS[1]
outputfilepath=ARGS[2]
conc_dim=parse(Int,ARGS[3])

f = open(inputfilepath)
firstit=true
dataX=[]
dataY=[]
dataZ=[]
dataV=[]
for line in eachline(f)
    r = strip(line, ['\n'])
    print(r,"\n")
    dataXi = h5read(r, "/X_1")
    dataYi = h5read(r, "/Y_1")
    dataZi = h5read(r, "/Z_1")
    dataVi = h5read(r, "/var_1")
    if (firstit)
        dataX=dataXi
        dataY=dataYi
        dataZ=dataZi
        dataV=dataVi
        firstit=false
    else
        data_shape = collect(size(dataVi))
        data_shape[conc_dim] = data_shape[conc_dim]-1
        dataX=cat(conc_dim,dataX, dataXi[1:data_shape[1],1:data_shape[2],1:data_shape[3]])
        dataY=cat(conc_dim,dataY, dataYi[1:data_shape[1],1:data_shape[2],1:data_shape[3]])
        dataZ=cat(conc_dim,dataZ, dataZi[1:data_shape[1],1:data_shape[2],1:data_shape[3]])
        dataV=cat(conc_dim,dataV, dataVi[1:data_shape[1],1:data_shape[2],1:data_shape[3]])
    end
end
h5write(outputfilepath, "/X_1", dataX)
h5write(outputfilepath, "/Y_1", dataY)
h5write(outputfilepath, "/Z_1", dataZ)
h5write(outputfilepath, "/var_1", dataV)

print("total size = ",size(dataV),"\n")
print("change this in the accompanying xmf file...\n")
