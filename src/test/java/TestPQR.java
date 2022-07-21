import org.rcsb.biozernike.structure.*;
import org.rcsb.biozernike.volume.MapFileType;
import org.rcsb.biozernike.volume.VolumeIO;

import java.io.File;

import org.rcsb.biozernike.volume.Volume;



//input /home/org305/biozernike/DATASET/O94760/01_out/pockets/pocket1_vert.pqr

public class TestPQR {
    public static void main(String[] args) {
        double gridWidth = 0.1;
        StructurePQR structure = new StructurePQR(args[0], gridWidth);
        try{ 
            int[] dims = structure.calcBoundingBox(5);
            System.out.println("Dims: " + dims[0]+ " " + dims[1]+ " " + dims[2] + "\n");
            double[] voxels = new double[dims[0]*dims[1]*dims[2]];
            Volume volume = new Volume();
            volume.createFromData(dims, voxels,gridWidth);
            volume.resetVoxels();
            structure.fillVoxels(volume, dims,5);
            volume.updateCenter();
            //VolumeIO.write(volume, new File("fromPQRtoMRC.mrc"), MapFileType.MRC); 

            double[] descriptors = structure.calcDescriptor(volume, Integer.parseInt(args[1]));
            structure.writeDescriptors(descriptors, args[2],Integer.parseInt(args[1]), args[3]);
        } catch (Exception e){
            e.printStackTrace();  
        }
    }
    
}
