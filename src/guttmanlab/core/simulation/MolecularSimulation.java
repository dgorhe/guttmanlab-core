package guttmanlab.core.simulation;

import org.apache.commons.math3.ml.distance.EuclideanDistance;


	//for information about diffusion coefficients, see: http://book.bionumbers.org/what-are-the-time-scales-for-diffusion-in-cells/
	//time step = distance / (6 * diffusion coefficient )
	//solve for distance (e.g., mean distance traveled in time step t)
	//distance = sqrt( 6 * diffusion coefficient )
	//the above equation will be used to calculate the mean distance traveled in one time step given a diffusion coefficient 

public class MolecularSimulation {

	//define input simulation parameters
		double outside_diffusion_coefficient = 0.05;
		int moleculeCopyNumber = 10000;
		int time_steps = 1000;
		double radius_of_nucleus = 2.5;
		double radius_of_inner_compartment = 0.7; // Xist occupies ~1/46th of the volume of the nucleus
		//double radius_of_inner_compartment = 0.35; // Xist occupies ~1/46th of the volume of the nucleus
		double[] inner_compartment_xyz = {0,0,0};
		
	//run one simulation given a ratio of outside to inside diffusion coefficients
		private double runSimulation(double insideToOutsideRatio){
			//parameters derived from inputs
			double outside_distance_per_step = Math.sqrt(6 * outside_diffusion_coefficient);
			double inside_diffusion_coefficient = outside_diffusion_coefficient / insideToOutsideRatio;
			double inside_distance_per_step = Math.sqrt(6 * inside_diffusion_coefficient);
			double[][] moleculePositions = new double[this.moleculeCopyNumber][];
			
			for (int i=0; i<this.moleculeCopyNumber; i++){
				moleculePositions[i]=randomPointInSphere(radius_of_nucleus);
			}
				
			
			int finalMoleculesInside = 0;
			
			for(int step=0; step<time_steps; step++){
				int numMoleculesInside = 0;
				
				//update position of each molecule
				for(int molecule=0; molecule<moleculeCopyNumber; molecule++){
					double[] moleculeXYZ = moleculePositions[molecule];
					double distancePerStep = outside_distance_per_step;
					
					// check if molecule is inside inner compartment
					double distance_from_inner_compartment = new EuclideanDistance().compute(moleculeXYZ, inner_compartment_xyz);
					if (distance_from_inner_compartment < radius_of_inner_compartment){
						//if molecule is inside inner compartment, use inside_distance_per_step
						distancePerStep = inside_distance_per_step;
						numMoleculesInside += 1;
					}
					
					//update molecule coordinates
					double phi = Math.random()*(Math.PI*2);
					double costheta = (Math.random()*2.0)-1.0;
					double theta = Math.acos( costheta );
					double x = Math.sin( theta ) * Math.cos( phi ) * distancePerStep;
					double y = Math.sin( theta ) * Math.sin( phi ) * distancePerStep;
					double z = Math.cos( theta ) * distancePerStep;
					double[] updated_position = {moleculeXYZ[0] + x, moleculeXYZ[1] + y, moleculeXYZ[2] + z};
				
					//check if updated position is inside nucleus
					if (Math.sqrt(Math.pow(updated_position[0], 2) + Math.pow(updated_position[1], 2) + Math.pow(updated_position[2], 2)) < radius_of_nucleus){
						moleculePositions[molecule] = updated_position;
					}			
				}
				finalMoleculesInside = numMoleculesInside;
			}
					
			return (double)finalMoleculesInside / (double)moleculeCopyNumber;
		}
	
	//source: https://karthikkaranth.me/blog/generating-random-points-in-a-sphere/
		private double[] randomPointInSphere(double radius){
			double u = Math.random();
			double v = Math.random();
			double theta = u * 2.0 * Math.PI;
			double phi = Math.acos(2.0 * v - 1.0);
			double r = Math.pow(Math.random(), (1.0 / 3)) * radius;
			double sintheta = Math.sin( theta );
			double costheta = Math.cos( theta );
			double sinphi = Math.sin( phi );
			double cosphi = Math.cos( phi );
			double x = r * sinphi * costheta;
			double y = r * sinphi * sintheta;
			double z = r * cosphi;
			double[] rtrn={x,y,z};
			return rtrn;
		}
	
		public static void main(String[] args){
			/*double[] ratios = {1,2,4,8,16,32,64,128,256,512,1024};
			for(int i=0; i<ratios.length; i++){
				double fraction=new MolecularSimulation().runSimulation(ratios[i]);
				System.err.println(ratios[i]+" "+fraction);
			}*/
			
			//double ratioOfHX=10000.0;
			double KdH=Math.pow(10, -8); //SAHA
			//double KdH=Math.pow(10, -6); //MS-275
			double KdX= Math.pow(10, -9);
			
			double[] ratiosOfKd={1, 2, 10, 100, 1000, 10000, 100000, 1000000, Math.pow(10, 7), Math.pow(10, 8)};
			
			double fractionInitial=new MolecularSimulation().runSimulation(1);
			for(int i=0; i<ratiosOfKd.length; i++){
				double ratio=ratiosOfKd[i];
				double fraction=new MolecularSimulation().runSimulation(ratio);
				System.err.println(ratiosOfKd[i]+" "+fraction+" "+(fraction/fractionInitial));
			}
			
		}
		
}
