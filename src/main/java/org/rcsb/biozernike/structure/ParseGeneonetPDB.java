package org.rcsb.biozernike.structure;

import org.rcsb.biozernike.volume.Volume;
import org.biojava.nbio.structure.Atom;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.List;
import org.rcsb.biozernike.descriptor.Descriptor;
import org.rcsb.biozernike.descriptor.DescriptorConfig;
import org.rcsb.biozernike.descriptor.DescriptorMode;

public class ParseGeneonetPDB {

    public int[] calcBoundingBox(Atom[] atoms, double[] gridWidths, int padding){
        
        int [] dims = new int[3];
        double[] min = getMin(atoms);
        double[] max = getMax(atoms);

        for(int i = 0; i<3; i++){
            dims[i] = (int)(((max[i]-min[i])/gridWidths[i])+2*padding);
        }
        return dims;

    }

    public double[] getMin(Atom[] atoms){
        double min_x = atoms[0].getX();
        double min_y = atoms[0].getY();
        double min_z = atoms[0].getZ();

        for(int i = 0; i<atoms.length; i++){
            if(atoms[i].getX()<min_x){
                min_x = atoms[i].getX();
            }
            if(atoms[i].getY()<min_y){
                min_y = atoms[i].getY();
            }
            if(atoms[i].getZ()<min_z){
                min_z = atoms[i].getZ();
            }
        }
        return new double[] {min_x, min_y, min_z};
    }


    public double[] getMax(Atom[] atoms){
        double max_x = atoms[0].getX();
        double max_y = atoms[0].getY();
        double max_z = atoms[0].getZ();

        for(int i = 0; i<atoms.length; i++){
            if(atoms[i].getX()>max_x){
                max_x = atoms[i].getX();
            }
            if(atoms[i].getY()>max_y){
                max_y = atoms[i].getY();
            }
            if(atoms[i].getZ()>max_z){
                max_z = atoms[i].getZ();
            }    
            
        }
        return new double[]{max_x, max_y, max_z};
    }

    public double[] calcGridWidths(Atom[] atoms){
        double diffx, diffy, diffz;
        double [] widths = new double[] {10,10,10};
        for(int i = 0; i<atoms.length-1; i++){
            diffx = Math.abs(atoms[i+1].getX()-atoms[i].getX());
            diffy = Math.abs(atoms[i+1].getY()-atoms[i].getY());
            diffz = Math.abs(atoms[i+1].getZ()-atoms[i].getZ());
            if(diffx < widths[0] && diffx != 0){
                widths[0]=diffx;
            }
            if(diffy < widths[1] && diffy != 0){
                widths[1]=diffy;
            }
            if(diffz < widths[2] && diffz != 0){
                widths[2]=diffz;
            }
        }
        return widths;
    }

    public void fillVoxels(Atom[] atoms, double[] gridWidths, Volume volume, int[] dims, int padding){
        double[] min = getMin(atoms);
        int idx, idy, idz;
            for(Atom atom:atoms){
                idx = (int)Math.round((atom.getX()-min[0])/gridWidths[0])+padding;
                idy = (int)Math.round((atom.getY()-min[1])/gridWidths[1])+padding;
                idz = (int)Math.round((atom.getZ()-min[2])/gridWidths[2])+padding;
                volume.setValue(idx, idy, idz, 1.0);
            }
    }


    public void writeDescriptors(double[] descriptors, String nameFile){
        FileWriter fileWriter;
        try {
            fileWriter = new FileWriter(nameFile);
            for(int j = 0; j<descriptors.length; j++){
                fileWriter.write(j + " " + descriptors[j]+ "\n");
            }
            fileWriter.close();
        } catch (IOException e) {
             e.printStackTrace();
        }
    }


    public double [] calcDescriptor(Volume volume, int order){
		DescriptorConfig config;
        config = new DescriptorConfig(order, new int[] {0});
        Descriptor ssd = new Descriptor(volume,config);
        List<List<Double>> desc = ssd.getMomentInvariantsRaw();
        double[] di = desc.get(0).stream().mapToDouble(d -> d).toArray();
        return di;
        
    }

}
