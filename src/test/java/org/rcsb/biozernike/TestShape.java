package org.rcsb.biozernike;

import org.rcsb.biozernike.volume.Volume;
import org.jheaps.DoubleEndedAddressableHeap;
import org.junit.Test;
import org.rcsb.biozernike.volume.MapFileType;
import org.rcsb.biozernike.volume.VolumeIO;
import org.rcsb.biozernike.zernike.ZernikeMoments;

import static org.junit.Assert.assertArrayEquals;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.EnumSet;
import java.util.List;
import org.rcsb.biozernike.complex.Complex;
import org.rcsb.biozernike.descriptor.Descriptor;
import org.rcsb.biozernike.descriptor.DescriptorConfig;
import org.rcsb.biozernike.descriptor.DescriptorMode;

public class TestShape {
    
    @Test
    public void testCube(){
        int [] dims = new int[] {60,60,60};
        double[] voxels = new double[dims[0]*dims[1]*dims[2]];
        Volume volume = new Volume();
        volume.createFromData(dims, voxels,1.0);
        volume.resetVoxels();
        createCube(dims, 3.0, volume);
        System.out.println("Moments for Cube\n");
        calculateMomentsAndWriteVolume(volume, "shapes/cube.mrc");
        System.out.println("Descriptors for Cube\n");
        calcDescriptor(volume);
    }

    @Test
    public void testCube1(){
        int [] dims = new int[] {100,100,100};
        double[] voxels = new double[dims[0]*dims[1]*dims[2]];
        Volume volume = new Volume();
        volume.createFromData(dims, voxels,1.0);
        volume.resetVoxels();
        createCube(dims, 26.0, volume);
        System.out.println("Moments for cube1\n");
        calculateMomentsAndWriteVolume(volume, "shapes/cube1.mrc");
        System.out.println("Descriptors for cube1\n");
        calcDescriptor(volume);
        
    }

    public double[] AllCubes(int l){
        int [] dims = new int[] {100,100,100};
        double[] voxels = new double[dims[0]*dims[1]*dims[2]];
        Volume volume = new Volume();
        volume.createFromData(dims, voxels,1.0);
        volume.resetVoxels();
        createCube(dims, l, volume);
        volume.updateCenter();
        return calcDescriptor(volume);
    }


    public void Loop(){
        FileWriter fileWriter;
        try {
            fileWriter = new FileWriter("cubes.txt");
            for(int i = 2; i<50; i+=2){
                double[] descriptors = AllCubes(i);
                for(int j = 0; j<descriptors.length; j++){
                    fileWriter.write(i + " " + j + " " + descriptors[j]+ "\n");
                }
                System.out.println("Ho scritto il cubo " + i + "\n");
            }
            fileWriter.close();
        } catch (IOException e) {
             e.printStackTrace();
        }
    }


    public double[] AllCylinders(double r, double h){
        int [] dims = new int[] {100,100,100};
        double[] voxels = new double[dims[0]*dims[1]*dims[2]];
        Volume volume = new Volume();
        volume.createFromData(dims, voxels,1.0);
        volume.resetVoxels();
        createCylinder(dims,r,h, volume);
        volume.updateCenter();
        return calcDescriptor(volume);
    }

    public void LoopCylinder(){
        FileWriter fileWriter;
        try {
            fileWriter = new FileWriter("cylinders.txt");
            for(int i = 1; i<60; i+=2){
                for(int k = 1; k<50; k+=4){
                    double[] descriptors = AllCylinders(i, (i+3)/2);
                    for(int j = 0; j<descriptors.length; j++){
                        fileWriter.write(i + " " + k + " " + j +  " " + descriptors[j]+ "\n");
                    }
                }
                System.out.println("Ho scritto il cilindro " + i + "\n");
            }
            fileWriter.close(); 
        } catch (IOException e) {
             e.printStackTrace();
        }
    }


    public Volume rotateShapeXZ(double angle, Volume volume, int n){ 
        double[][] rot = {{Math.cos(angle),0, -Math.sin(angle)},
                        {0,1,0},
                        {Math.sin(angle),0, Math.cos(angle)}};
        int [] dims = new int[] {n,n,n};
        double[] voxels = new double[n*n*n];
        Volume vrot = new Volume();
        vrot.createFromData(dims, voxels,1.0);
        vrot.resetVoxels();
            
        double xp,yp,zp;
        double xm, ym, zm;
        //[x', y', z']=rot*[x,y,z] + t
        for(int z = 0; z<n; z++){
            for(int y = 0; y<n; y++){
                for(int x = 0; x<n; x++){
                    xm = x-n/2;
                    ym = y-n/2;
                    zm = z-n/2;
                    if(volume.getValue(x, y, z)==1){
                        xp = (xm*rot[0][0]+ym*rot[0][1]+zm*rot[0][2]) + n/2;
                        yp = (xm*rot[1][0]+ym*rot[1][1]+zm*rot[1][2]) + n/2;
                        zp = (xm*rot[2][0]+ym*rot[2][1]+zm*rot[2][2]) + n/2;
                        vrot.setValue((int)xp, (int)yp, (int)zp, 1);  
                    }
                }
            }
        }
        return vrot;
    }

    public void LoopRotate(){
        FileWriter fileWriter;
        double[] angles = new double[]{0,6,12,18,24,30,36,42,48,54,60,66,72,78,84,90};
        int [] dims = new int[] {150,150,150};
        double[] voxels = new double[dims[0]*dims[1]*dims[2]];
        Volume volume = new Volume();
        volume.createFromData(dims, voxels,1.0);
        volume.resetVoxels();
        createCylinder(dims,27.0,40.0, volume);
        volume.updateCenter();

        try {
            fileWriter = new FileWriter("cylindersRotate.txt");
            for(int i = 0; i<angles.length; i++){
                Volume vrot = rotateShapeXZ((angles[i]*Math.PI)/180,volume, dims[0]);
                calculateMomentsAndWriteVolume(vrot, "shapes/vrot" + angles[i]+".mrc");
                double[] descriptors = calcDescriptor(vrot);
                for(int j = 0; j<descriptors.length; j++){
                    fileWriter.write(angles[i] + " " + j +  " " + descriptors[j]+ "\n");
                }
                System.out.println("Ho scritto il cilindro ruotato con angolo " + angles[i] + "\n");
            }
            fileWriter.close(); 
        } catch (IOException e) {
             e.printStackTrace();
        }
    }

    public double[] AllSpheres(double r){
        int [] dims = new int[] {100,100,100};
        double[] voxels = new double[dims[0]*dims[1]*dims[2]];
        Volume volume = new Volume();
        volume.createFromData(dims, voxels,1.0);
        volume.resetVoxels();
        createSphere(dims, r, volume);
        volume.updateCenter();
        return calcDescriptor(volume);
    }

    public void LoopSphere(){
        FileWriter fileWriter;
        try {
            fileWriter = new FileWriter("spheres.txt");
            for(int i = 1; i<60; i+=2){
                double[] descriptors = AllSpheres(i);
                for(int j = 0; j<descriptors.length; j++){
                    fileWriter.write(i + " " + j + " " + descriptors[j]+ "\n");
                }
                System.out.println("Ho scritto la sfera " + i + "\n");
            }
            fileWriter.close(); 
        } catch (IOException e) {
             e.printStackTrace();
        }
    }


    @Test
    public void testCube2(){
        int [] dims = new int[] {100,100,100};
        double[] voxels = new double[dims[0]*dims[1]*dims[2]];
        Volume volume = new Volume();
        volume.createFromData(dims, voxels,1.0);
        volume.resetVoxels();
        createCube(dims, 38.0, volume);
        System.out.println("Moments for cube2\n");
        calculateMomentsAndWriteVolume(volume, "shapes/cube2.mrc");
        System.out.println("Descriptors for cube2\n");
        calcDescriptor(volume);
    }



    @Test
    public void testCylinder(){ 
        int [] dims = new int[] {80,80,80};
        double[] voxels = new double[dims[0]*dims[1]*dims[2]];
        Volume volume = new Volume();
        volume.createFromData(dims, voxels,1.0);
        volume.resetVoxels();
        createCylinder(dims, 3.0,9.0, volume);
        System.out.println("Moments for cylinder\n");
        calculateMomentsAndWriteVolume(volume, "shapes/cylinder.mrc");
        System.out.println("Descriptors for cylinder\n");
        calcDescriptor(volume);
    }

    @Test
    public void testCylinder1(){ 
        int [] dims = new int[] {150,150,150};
        double[] voxels = new double[dims[0]*dims[1]*dims[2]];
        Volume volume = new Volume();
        volume.createFromData(dims, voxels,1.0);
        volume.resetVoxels();
        createCylinder(dims, 15.0,37.0, volume);
        System.out.println("Moments for cylinder1\n");
        List<List<List<Complex>>> momentcylinder = calculateMomentsAndWriteVolume(volume, "shapes/cylinder1.mrc");
        System.out.println("Descriptors for cylinder1\n");
        double[] descriptorCylinder = calcDescriptor(volume);

       //translate

        Volume vtran = translate(dims[0], volume);
        System.out.println("Moments for cylinder (translate)");
        List<List<List<Complex>>> momentTran = calculateMomentsAndWriteVolume(vtran, "shapes/cilindroTraslato1.mrc");
        System.out.println("Descriptors for cylinder (translate)");
        double [] descriptorTranslate = calcDescriptor(vtran);

        //assertArrayEquals(descriptorCylinder, descriptorTranslate,1e-9);

      
        
    }

    @Test
    public void testSphere(){ 
        int [] dims = new int[] {100,100,100};
        double[] voxels = new double[dims[0]*dims[1]*dims[2]];
        Volume volume = new Volume();
        volume.createFromData(dims, voxels,1.0);
        volume.resetVoxels();
        createSphere(dims, 8.0, volume);
        System.out.println("Moments for sphere\n");
        calculateMomentsAndWriteVolume(volume, "shapes/sphere.mrc");
        System.out.println("Descriptors for sphere\n");
        calcDescriptor(volume);
    }

    public void testSphereNegative(){ 
        int [] dims = new int[] {100,100,100};
        double[] voxels = new double[dims[0]*dims[1]*dims[2]];
        Volume volume = new Volume();
        volume.createFromData(dims, voxels,1.0);
        volume.resetVoxels();
        createNegativeSphere(dims, 30.0, volume);
        System.out.println("Moments for sphere negative\n");
        calculateMomentsAndWriteVolume(volume, "shapes/sphereNegative.mrc");
        System.out.println("Descriptors for sphere negative\n");
        calcDescriptor(volume);
    }


    public void testSphere1(){ 
        int [] dims = new int[] {100,100,100};
        double[] voxels = new double[dims[0]*dims[1]*dims[2]];
        Volume volume = new Volume();
        volume.createFromData(dims, voxels,1.0);
        //volume.setAllVoxels(-1.0);
        volume.resetVoxels();
        createSphere(dims, 30.0, volume);
        System.out.println("Moments for sphere1\n");
        calculateMomentsAndWriteVolume(volume, "shapes/sphere1.mrc");
        System.out.println("Descriptors for sphere1\n");
        calcDescriptor(volume);
    }

    public void testSphere2(){
        int [] dims = new int[] {100,100,100};
        double[] voxels = new double[dims[0]*dims[1]*dims[2]];
        Volume volume = new Volume();
        volume.createFromData(dims, voxels,1.0);
        volume.resetVoxels();
        createSphere(dims, 12.0, volume);
        System.out.println("Moments fot sphere2\n");
        calculateMomentsAndWriteVolume(volume, "shapes/sphere2.mrc");
        System.out.println("Descriptors for sphere2\n");
        calcDescriptor(volume);
    }

    public  List<List<List<Complex>>> calculateMomentsAndWriteVolume(Volume volume, String nameFile)  {
		volume.updateCenter(); 
        int i = 0;
        try{
            VolumeIO.write(volume, new File(nameFile), MapFileType.MRC); 
        }catch(IOException e){
            System.out.println("Error");
        }
        ZernikeMoments z = new ZernikeMoments(volume, 10); 
        List<List<List<Complex>>> originalMoments = z.getOriginalMoments();
        for(List<List<Complex>> moment :originalMoments){
            System.out.println("order " + i + " :" + moment);
            i++;
        }
        return originalMoments;
    }

    public double [] calcDescriptor(Volume volume){
        
        EnumSet<DescriptorMode> mode = EnumSet.of(DescriptorMode.CALCULATE_RAW);
        int i;
		DescriptorConfig config;
        //config = new DescriptorConfig(DescriptorTest.class.getResourceAsStream("/descriptor.properties"), mode);
        config = new DescriptorConfig(20, new int[] {0});
        Descriptor ssd = new Descriptor(volume,config);
        List<List<Double>> desc = ssd.getMomentInvariantsRaw();
        double[] di = desc.get(0).stream().mapToDouble(d -> d).toArray();
        return di;
        
    }


    public void createCube(int[] dims, double l, Volume volume){
        for(int z = 0; z < dims[2]; z++){
            for(int y = 0; y < dims[1]; y++){
                for(int x = 0; x < dims[0]; x++){
                    if(( x >= dims[0]/2 - l/2 && x <= dims[0]/2 + l/2) && ( y >= dims[1]/2 - l/2 && y <= dims[1]/2 + l/2) && ( z >= dims[2]/2 - l/2 && z <= dims[2]/2 + l/2)){
                        volume.setValue(x, y, z, 1.0);
                    }
                }

            }
        }

    }

    public void createSphere(int[] dims, double r, Volume volume){
        for(int z = 0; z < dims[2]; z++){
            for(int y = 0; y < dims[1]; y++){
                for(int x = 0; x < dims[0]; x++){
                    if( Math.sqrt(Math.pow(x-dims[0]/2,2)+Math.pow(y-dims[1]/2,2)+ Math.pow(z-dims[2]/2, 2)) <= r ){ 
                        volume.setValue(x, y, z, 1.0);
                    }
                }

            }
        }
    }

    public void createNegativeSphere(int[] dims, double r, Volume volume){
        for(int z = 0; z < dims[2]; z++){
            for(int y = 0; y < dims[1]; y++){
                for(int x = 0; x < dims[0]; x++){
                    if( Math.sqrt(Math.pow(x-dims[0]/2,2)+Math.pow(y-dims[1]/2,2)+ Math.pow(z-dims[2]/2, 2)) <= r ){ 
                        volume.setValue(x, y, z, -2.0);
                    }
                }

            }
        }
    }

    public void createCylinder(int[] dims, double r, double h, Volume volume){
        for(int z = 0; z < dims[2]; z++){
            for(int y = 0; y < dims[1]; y++){
                for(int x = 0; x < dims[0]; x++){
                    if( z>= (dims[2]/2 - h/2) && z <= (dims[2]/2) + h/2){ 
                        if(Math.sqrt(Math.pow(x-dims[0]/2,2) + Math.pow(y-dims[1]/2,2)) <= r){
                            volume.setValue(x, y, z, 1.0);
                        }
                       
                    }
                }

            }
        }
    }


    public Volume translate(int n, Volume volume){ 
        int t[] = {10,10,10}; //vettore di traslazione
        int [] dims = new int[] {n,n,n};
        double[] voxels = new double[n*n*n];
        Volume vtran = new Volume();
        vtran.createFromData(dims, voxels,1.0);
        vtran.resetVoxels();
        
        int xp,yp,zp;
        //[x', y', z']=[x,y,z]+t
        for(int z = 0; z<n-t[2]; z++){
            for(int y = 0; y<n-t[1]; y++){
                for(int x = 0; x<n-t[0]; x++){
                    if(volume.getValue(x, y, z)==1){
                        xp = x+t[0];
                        yp = y+t[1];
                        zp = z+t[2];
                        vtran.setValue(xp, yp, zp, volume.getValue(x, y, z));  
                    }
                }
            }
        }
        return vtran;
    }


    public static void main(String[] args) {
        TestShape t = new TestShape();
        /*t.testCube();
        t.testCube1();
        t.testCube2();
        t.testSphere();
        t.testSphere1();
        t.testSphere2();*/
        //t.testSphereNegative();
        /*t.testCylinder();
        t.testCylinder1();*/
        /*t.Loop();
        t.LoopCylinder();
        t.LoopSphere();*/
        t.LoopRotate();
       
    }
}

