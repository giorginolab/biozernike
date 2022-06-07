import org.rcsb.biozernike.structure.*;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.io.PDBFileReader;
import org.rcsb.biozernike.volume.MapFileType;
import org.rcsb.biozernike.volume.VolumeIO;

import java.io.File;

import org.rcsb.biozernike.volume.Volume;
import org.biojava.nbio.structure.Atom;

//input /home/org305/biozernike/DATASET/O94760/01_h_pocket1.pdb

public class TestGENEONET {
    public static void main(String[] args) {
        PDBFileReader pdbreader = new PDBFileReader();
        ParseGeneonetPDB c = new ParseGeneonetPDB();
        try{
            Structure struc = pdbreader.getStructure(args[0]);  
            Atom[] atoms = StructureTools.getAllAtomArray(struc);
            double[] gridWidths = c.calcGridWidths(atoms);
            int[] dims = c.calcBoundingBox(atoms, gridWidths,5);
            System.out.println("Dims: " + dims[0]+ " " + dims[1]+ " " + dims[2] + "\n");

            double[] voxels = new double[dims[0]*dims[1]*dims[2]];
            Volume volume = new Volume();
            volume.createFromData(dims, voxels,(double)(gridWidths[0]+gridWidths[1]+gridWidths[2])/3);
            volume.resetVoxels();
            c.fillVoxels(atoms, gridWidths, volume, dims,5);
            volume.updateCenter();

            VolumeIO.write(volume, new File("fromPDBtoMRC.mrc"), MapFileType.MRC); 
            double[] descriptors = c.calcDescriptor(volume, Integer.parseInt(args[1]));
            c.writeDescriptors(descriptors, args[2], Integer.parseInt(args[1]));
        } catch (Exception e){
            e.printStackTrace();  
        }
    }
}
