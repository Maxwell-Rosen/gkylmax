import h5py
import numpy as np
import matplotlib.pyplot as plt

hdf5_file = 'cql3d_f_RZ.h5'  # Replace with your HDF5 file path
output_dir = '../binary_files'  # Output directory for binary files

def convert_hdf5_to_binary(hdf5_file, output_dir):
    # Open the HDF5 file
    with h5py.File(hdf5_file, 'r') as f:
        # Iterate over attributes and save each one to a binary file
        for attr_name, attr_data in f.items():
            print(f"Converting attribute: {attr_name}")
            data = np.array(attr_data)
            dtype = data.dtype
            if dtype == '>f8':
                data = data.astype(np.float64)

            if attr_name == 'f_dist':
                # continue
                data_ion = data[0,:,:,:,:] * 1e12 # 0 is for ions
                data_electron = data[1,:,:,:,:] * 1e12 # 1 is for electrons

                output_file = f"{output_dir}/f_dist_ion.bin"
                # Write the data to binary file
                with open(output_file, 'wb') as binary_file:
                    binary_file.write(data_ion.tobytes())

                output_file = f"{output_dir}/f_dist_elc.bin"
                # Write the data to binary file
                with open(output_file, 'wb') as binary_file:
                    binary_file.write(data_electron.tobytes())
            elif attr_name == 'psiGrid':
                output_file = f"{output_dir}/psiGrid.bin"
                # Write the data to binary file
                data *= -1e-8
                with open(output_file, 'wb') as binary_file:
                    binary_file.write(data.tobytes())
            elif attr_name == 'v_norm':
                output_file = f"{output_dir}/v_norm.bin"
                # Write the data to binary file
                data *= 1e-2
                with open(output_file, 'wb') as binary_file:
                    binary_file.write(data.tobytes())
            elif attr_name == 'zGrid':
                output_file = f"{output_dir}/zGrid.bin"
                # Write the data to binary file
                data *= 1e-2
                with open(output_file, 'wb') as binary_file:
                    binary_file.write(data.tobytes())
            elif attr_name == 'vGrid':
                output_file = f"{output_dir}/vGrid.bin"
                # Write the data to binary file
                data *= 1e-2
                with open(output_file, 'wb') as binary_file:
                    binary_file.write(data.tobytes())
            elif attr_name == 'BdB0':
                output_file = f"{output_dir}/BdB0.bin"
                # Write the data to binary file
                with open(output_file, 'wb') as binary_file:
                    binary_file.write(data.tobytes())

                output_file = f"{output_dir}/BGrid.bin"
                ## Find the output B0 attribute
                B0 = f['B0'][()]
                tile_B0 = np.transpose(np.tile(B0, (data.shape[1], 1)))
                BGrid = np.multiply(data, tile_B0) * 1e-4
                # Write the data to binary file
                with open(output_file, 'wb') as binary_file:
                    binary_file.write(BGrid.tobytes())
            else:
                output_file = f"{output_dir}/{attr_name}.bin"
                # Write the data to binary file
                with open(output_file, 'wb') as binary_file:
                    binary_file.write(data.tobytes())
        





    # hf.create_dataset('v_norm',data=v_norm)    #[cm/s] * 1e-2
    # hf.create_dataset('sqPsiGrid',data=sqPsiGrid) #[normalized]
    # hf.create_dataset('psiLim',data=psiLim)    #[normalized]
    # hf.create_dataset('psiGrid',data=psiGrid)  #[Gauss cm^2] 
    # hf.create_dataset('zGrid',data=zGrid[0,:]) #[cm] * 1e-2
    # hf.create_dataset('uGrid',data=x)          #[normalized]
    # hf.create_dataset('vGrid',data=v0)         #[cm/s] * 1e-2
    # hf.create_dataset('theta',data=theta)      #[radians]
    # hf.create_dataset('charge',data=q)         #[e]
    # hf.create_dataset('mass',data=m)           #[g]
    # hf.create_dataset('BdB0',data=BdB0)        #[normalized]
    # hf.create_dataset('B0',data=B0)            #[Gauss]
    # hf.create_dataset('phi',data=phi)          #[kV]
    # hf.create_dataset('f_dist',data=f_out)     #[cm^-6 s^3] * 1e12

convert_hdf5_to_binary(hdf5_file, output_dir)
