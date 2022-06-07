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

    public StructurePQR(String s){
        x = new ArrayList<>();
        y = new ArrayList<>();
        z = new ArrayList<>();
        radius = new ArrayList<>();
        read(s);
    }

    public void read(String s){
        try{
            FileReader fr = new FileReader(s);
            BufferedReader br = new BufferedReader(fr);
            String line;
            String split = " ";

            while((line=br.readLine())!=null){
                String[] parts = line.split(split);
                if(parts[0].equals("ATOM")){
                    if(parts[5].isEmpty()){
                        x.add(Double.parseDouble(parts[22]));
                        y.add(Double.parseDouble(parts[25]));
                        z.add(Double.parseDouble(parts[27]));
                        radius.add(Double.parseDouble(parts[36]));
                    }else{
                        x.add(Double.parseDouble(parts[21]));
                        y.add(Double.parseDouble(parts[24]));
                        z.add(Double.parseDouble(parts[26]));
                        radius.add(Double.parseDouble(parts[35]));
                    }
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
            dims[i] = (int)((max[i]-min[i])+2*padding);
        }
        return dims;
    }

    public void fillVoxels(Volume volume, int[] dims, int padding){
        double[] min = getMin();
        double idx, idy, idz;
        for(int i = 0; i<dims[0]; i++){
            for(int j = 0; j<dims[1]; j++){
                for(int k = 0; k<dims[2]; k++){
                    for(int h = 0; h< x.size(); h++){
                        idx = i-padding + min[0];
                        idy = j-padding + min[1];
                        idz = k-padding + min[2];
                        if( Math.sqrt(Math.pow(idx-x.get(h),2)+Math.pow(idy-y.get(h),2)+ Math.pow(idz-z.get(h), 2)) <= radius.get(h) ){ 
                            volume.setValue(i,j,k, 1.0);
                        }
                    }
                }
            }
        }
    }

    public void writeDescriptors(double[] descriptors, String nameFile, int order){
        FileWriter fileWriter;
        int j = 0;

        try {
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


    
    
    



