package org.rcsb.biozernike;

import org.rcsb.biozernike.volume.Volume;
import org.apache.commons.lang.ArrayUtils;
import org.junit.Test;
import org.rcsb.biozernike.volume.MapFileType;
import org.rcsb.biozernike.volume.VolumeIO;
import org.rcsb.biozernike.zernike.ZernikeMoments;

import static org.junit.Assert.assertArrayEquals;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.activation.FileTypeMap;
import javax.vecmath.Matrix4d;

import org.rcsb.biozernike.complex.Complex;

public class TestShape {
    
    @Test
    public void testCube(){
        int [] dims = new int[] {60,60,60};
        double[] voxels = new double[dims[0]*dims[1]*dims[2]];
        Volume volume = new Volume();
        volume.createFromData(dims, voxels,1.0);
        volume.resetVoxels();
        createCube(dims, 3.0, volume);
        calculateMomentsAndWriteVolume(volume, "shapes/cube.mrc");
    }

    @Test
    public void testCube1(){
        int [] dims = new int[] {60,60,60};
        double[] voxels = new double[dims[0]*dims[1]*dims[2]];
        Volume volume = new Volume();
        volume.createFromData(dims, voxels,1.0);
        volume.resetVoxels();
        createCube(dims, 21.0, volume);
        calculateMomentsAndWriteVolume(volume, "shapes/cube1.mrc");
        
    }


    @Test
    public void testCube2(){
        int [] dims = new int[] {60,60,60};
        double[] voxels = new double[dims[0]*dims[1]*dims[2]];
        Volume volume = new Volume();
        volume.createFromData(dims, voxels,1.0);
        volume.resetVoxels();
        createCube(dims, 38.0, volume);
        calculateMomentsAndWriteVolume(volume, "shapes/cube2.mrc");
    }



    @Test
    public void testCylinder(){ 
        int [] dims = new int[] {80,80,80};
        double[] voxels = new double[dims[0]*dims[1]*dims[2]];
        Volume volume = new Volume();
        volume.createFromData(dims, voxels,1.0);
        volume.resetVoxels();
        createCylinder(dims, 3.0,9.0, volume);
    }

    @Test
    public void testCylinder1(){ 
        int [] dims = new int[] {80,80,80};
        double[] voxels = new double[dims[0]*dims[1]*dims[2]];
        Volume volume = new Volume();
        volume.createFromData(dims, voxels,1.0);
        volume.resetVoxels();
        createCylinder(dims, 12.0,51.0, volume);
        List<List<List<Complex>>> momentcylinder = calculateMomentsAndWriteVolume(volume, "shapes/cylinder.mrc");

        //rotate
        Volume vrot = rotate(dims[0], volume);
        System.out.println("Descriptors for cylinder (rotate)");
        List<List<List<Complex>>> momentRot = calculateMomentsAndWriteVolume(vrot, "shapes/cilindroRuotato.mrc");

        double[] momentcylinderArr = ArrayUtils.toPrimitive(
				ZernikeMoments.flattenMomentsDouble(
						momentcylinder).
						toArray(new Double[0]));

        double[] momentRotArr = ArrayUtils.toPrimitive(
                ZernikeMoments.flattenMomentsDouble(
                        momentRot).
                        toArray(new Double[0]));
       
        
       //assertArrayEquals(momentcylinderArr,momentRotArr,1e-10);

       //translate

       Volume vtran = translate(dims[0], volume);
        System.out.println("Descriptors for cylinder (translate)");
        List<List<List<Complex>>> momentTran = calculateMomentsAndWriteVolume(vtran, "shapes/cilindroTraslato.mrc");


        double[] momentTranArr = ArrayUtils.toPrimitive(
                ZernikeMoments.flattenMomentsDouble(
                        momentTran).
                        toArray(new Double[0]));
        
        assertArrayEquals(momentcylinderArr,momentTranArr,1e-10);
        
    }

    @Test
    public void testSphere(){ 
        int [] dims = new int[] {100,100,100};
        double[] voxels = new double[dims[0]*dims[1]*dims[2]];
        Volume volume = new Volume();
        volume.createFromData(dims, voxels,1.0);
        volume.resetVoxels();
        createSphere(dims, 8.0, volume);
        calculateMomentsAndWriteVolume(volume, "shapes/sphere.mrc");
    }


    public void testSphere1(){ 
        int [] dims = new int[] {100,100,100};
        double[] voxels = new double[dims[0]*dims[1]*dims[2]];
        Volume volume = new Volume();
        volume.createFromData(dims, voxels,1.0);
        volume.resetVoxels();
        createSphere(dims, 3.0, volume);
        calculateMomentsAndWriteVolume(volume, "shapes/sphere1.mrc");
    }

    public void testSphere2(){
        int [] dims = new int[] {100,100,100};
        double[] voxels = new double[dims[0]*dims[1]*dims[2]];
        Volume volume = new Volume();
        volume.createFromData(dims, voxels,1.0);
        volume.resetVoxels();
        createSphere(dims, 12.0, volume);
        calculateMomentsAndWriteVolume(volume, "shapes/sphere2.mrc");
    }

    public  List<List<List<Complex>>> calculateMomentsAndWriteVolume(Volume volume, String nameFile)  {
		volume.updateCenter(); 
        try{
            VolumeIO.write(volume, new File(nameFile), MapFileType.MRC); 
        }catch(IOException e){
            System.out.println("Error");
        }
        ZernikeMoments z = new ZernikeMoments(volume, 3); 
        List<List<List<Complex>>> originalMoments = z.getOriginalMoments();
        System.out.println(ZernikeMoments.flattenMomentsDouble(originalMoments));
        return originalMoments;
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

    public Volume rotate(int n, Volume volume){ 
        int[][] rot = {{1,0,0},{0,0,1},{0,1,0}};
        int [] dims = new int[] {n,n,n};
        double[] voxels = new double[n*n*n];
        Volume vrot = new Volume();
        vrot.createFromData(dims, voxels,1.0);
        vrot.resetVoxels();
        
        int xp,yp,zp;
        //[x', y', z']=rot*[x,y,z]
        for(int z = 0; z<n; z++){
            for(int y = 0; y<n; y++){
                for(int x = 0; x<n; x++){
                        if(volume.getValue(x, y, z)==1){
                            xp = x*rot[0][0]+y*rot[0][1]+z*rot[0][2];
                            yp = x*rot[1][0]+y*rot[1][1]+z*rot[1][2];
                            zp = x*rot[2][0]+y*rot[2][1]+z*rot[2][2];
                            vrot.setValue(xp, yp, zp, volume.getValue(x, y, z));  
                        }
                }
            }
        }
        return vrot;
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
        /*System.out.println("Descriptors for cube:");
        t.testCube();
        System.out.println("Descriptors for cube1:");
        t.testCube1();
        System.out.println("Descriptors for cube 2:");
        t.testCube2();
        System.out.println("Descriptors for sphere");
        t.testSphere();
        System.out.println("Descriptors for sphere1");
        t.testSphere1();
        System.out.println("Descriptors for sphere2");
        t.testSphere2();
        System.out.println("Descriptors for cylinder");
        t.testCylinder();*/
        System.out.println("Descriptors for cylinder1");
        t.testCylinder1();
       
    }
}

