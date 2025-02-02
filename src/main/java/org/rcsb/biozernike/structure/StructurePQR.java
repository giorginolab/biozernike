package org.rcsb.biozernike.structure;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import org.rcsb.biozernike.descriptor.Descriptor;
import java.util.List;

import org.rcsb.biozernike.descriptor.DescriptorConfig;
import org.rcsb.biozernike.volume.Volume;

public class StructurePQR {
    private ArrayList<Double> x;
    private ArrayList<Double> y;
    private ArrayList<Double> z;
    private ArrayList<Double> radius;
    private double gridWidth;

    public StructurePQR(String s, double gridWidth){
        x = new ArrayList<>();
        y = new ArrayList<>();
        z = new ArrayList<>();
        radius = new ArrayList<>();
        this.gridWidth = gridWidth;
        read(s);
    }

    public void read(String s){
        try{
            FileReader fr = new FileReader(s);
            BufferedReader br = new BufferedReader(fr);
            String line;
            String split = "\\s+";

            while((line=br.readLine())!=null){
                String[] parts = line.split(split);
                
                if(parts[0].equals("ATOM")){
                    x.add(Double.parseDouble(parts[5]));
                    y.add(Double.parseDouble(parts[6]));
                    z.add(Double.parseDouble(parts[7]));
                    radius.add(Double.parseDouble(parts[9]));
                }
               
               
            }
        }catch(FileNotFoundException e){
            System.out.println("File not found");
            System.exit(0);
        } catch (NumberFormatException e) {
            System.out.println("Format error");
            System.exit(0);
            e.printStackTrace();
        } catch (IOException e) {
            System.out.println("File not found");
            System.exit(0);
        }

    }

    public double[] getMin(){
        double min_x = x.get(0)-radius.get(0);
        double min_y = y.get(0)-radius.get(0);
        double min_z = z.get(0)-radius.get(0);

        for(int i = 0; i<x.size(); i++){
            if((x.get(i)-radius.get(i))<min_x){
                min_x = x.get(i)-radius.get(i);
            }
            if((y.get(i)-radius.get(i))<min_y){
                min_y = y.get(i)-radius.get(i);
            }
            if((z.get(i)-radius.get(i))<min_z){
                min_z = z.get(i)-radius.get(i);
            }
        }
        return new double[] {min_x, min_y, min_z};
    }


    public double[] getMax(){
        double max_x = x.get(0)+radius.get(0);
        double max_y = y.get(0)+radius.get(0);
        double max_z = z.get(0)+radius.get(0);

        for(int i = 0; i<x.size(); i++){
            if((x.get(i)+radius.get(i))>max_x){
                max_x = x.get(i)+radius.get(i);
            }
            if((y.get(i)+radius.get(i))>max_y){
                max_y = y.get(i)+radius.get(i);
            }
            if((z.get(i)+radius.get(i))>max_z){
                max_z = z.get(i)+radius.get(i);
            }    
            
        }
        return new double[]{max_x, max_y, max_z};
    }

    public int[] calcBoundingBox(int padding){
        
        int [] dims = new int[3];
        double[] min = getMin();
        double[] max = getMax();

        for(int i = 0; i<3; i++){
            dims[i] = (int)(((max[i]-min[i])/gridWidth)+2*padding);
        }
        return dims;
    }

    public void fillVoxels(Volume volume, int[] dims, int padding){
        double[] min = getMin();
        double coordx, coordy, coordz;
        for(int i = 0; i<dims[0]; i++){ 
             for(int j = 0; j<dims[1]; j++){
                for(int k = 0; k<dims[2]; k++){
                    for(int h = 0; h< x.size(); h++){
                        coordx = (i-padding)*gridWidth + min[0];
                        coordy = (j-padding)*gridWidth + min[1];
                        coordz = (k-padding)*gridWidth + min[2];
                        if( Math.sqrt(Math.pow(coordx-x.get(h),2)+Math.pow(coordy-y.get(h),2)+ 
                            Math.pow(coordz-z.get(h), 2)) <= radius.get(h) ){ 
                            volume.setValue(i,j,k, 1.0);
                        }
                    }
                }
            }
        }
    }

    public void writeDescriptors(double[] descriptors, String nameFile, int order, String nameFile1){
        FileWriter fileWriter;
        FileWriter fileWriter1;
        int j = 0;

        try {
        	//write file with index, descriptor, n, l, z and one row file
            fileWriter = new FileWriter(nameFile);
            fileWriter.write("Index Descriptor n l z\n");
            for(int n = 0; n<=order; n++){
                for(int k = 0; k<=(n/2); k++){
                    int l = k*2+(int)((n%2!=0) ? 1 : 0);
                    fileWriter.write(j + " " 
                    + descriptors[j]+ " " 
                    + n + " " 
                    + l + " " +
                    "z"+n+","+l
                    +  "\n");
                    j++;
                }
            }
            fileWriter.close();

            
            //write one row file
            j = 0;
            fileWriter1 = new FileWriter(nameFile1);
            fileWriter1.write("Name ");
            for(int n = 0; n<order; n++){
            	for(int k = 0; k<=(n/2); k++){
            	 int l = k*2+(int)((n%2!=0) ? 1 : 0);
            	 fileWriter1.write("z"+n+","+l + " ");
            	}
            }
            fileWriter1.write("\n");
            fileWriter1.write(nameFile);
             for(int n = 0; n<=order; n++){
                for(int k = 0; k<=(n/2); k++){
                    int l = k*2+(int)((n%2!=0) ? 1 : 0);
                    fileWriter1.write(" " 
                    + descriptors[j]);
                    j++;
                }
            }
            fileWriter1.write("\n");
            fileWriter1.close();
            
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


    
    
    



