package guttmanlab.core.simulation;

import org.apache.commons.math3.ml.distance.EuclideanDistance;

public class SimulateBivalentMolecules {

	
	//define input simulation parameters
	double outside_diffusion_coefficient = 0.05;
	//outside_to_inside_ratio = 1
	int sharp_copy_number = 1000;
	int time_steps = 1000;
	double radius_of_nucleus = 2.5;
	double radius_of_inner_compartment = 0.7; // Xist occupies ~1/46th of the volume of the nucleus
	double[] inner_compartment_xyz = {0,0,0};
	
	public SimulateBivalentMolecules (){}
		
	
	//source: https://karthikkaranth.me/blog/generating-random-points-in-a-sphere/
	private double[] randomPointInSphere(double radius){
		double u = Math.random();
		double v = Math.random();
		double theta = u * 2.0 * Math.PI;
		double phi = Math.acos(2.0 * v - 1.0);
		double r = Math.pow(Math.random(), (1.0 / 3)) * radius; //TODO Confirm
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
	
	//run one simulation given a ratio of outside to inside diffusion coefficients
	private double runSimulation(double outside_to_inside_ratio) {
		double outside_distance_per_step = Math.sqrt(6 * outside_diffusion_coefficient);
		double inside_diffusion_coefficient = outside_diffusion_coefficient / outside_to_inside_ratio;
		double inside_distance_per_step = Math.sqrt(6 * inside_diffusion_coefficient);
		double[][] sharp_positions = new double[this.sharp_copy_number][];
		for (int i=0; i<this.sharp_copy_number; i++){
			sharp_positions[i]=randomPointInSphere(radius_of_nucleus);
		}
			
		
		int final_num_sharp_inside = 0;
		
		for(int step=0; step<time_steps; step++){
			int num_sharp_inside = 0;
			
			//update position of each sharp molecule
			for(int sharpMolecule=0; sharpMolecule<sharp_copy_number; sharpMolecule++){
				double[] sharp_xyz = sharp_positions[sharpMolecule];
				double sharp_distance_per_step = outside_distance_per_step;
				
				// check if sharp is inside inner compartment
				double distance_from_inner_compartment = new EuclideanDistance().compute(sharp_xyz, inner_compartment_xyz);
				if (distance_from_inner_compartment < radius_of_inner_compartment){
					//if sharp is inside inner compartment, use inside_distance_per_step
					sharp_distance_per_step = inside_distance_per_step;
					num_sharp_inside += 1;
				}
				
				//update sharp coordinates
				double phi = Math.random()*(Math.PI*2);
				double costheta = (Math.random()*2.0)-1.0;
				double theta = Math.acos( costheta );
				double x = Math.sin( theta ) * Math.cos( phi ) * sharp_distance_per_step;
				double y = Math.sin( theta ) * Math.sin( phi ) * sharp_distance_per_step;
				double z = Math.cos( theta ) * sharp_distance_per_step;
				double[] updated_position = {sharp_xyz[0] + x, sharp_xyz[1] + y, sharp_xyz[2] + z};
			
				//check if updated position is inside nucleus
				if (Math.sqrt(Math.pow(updated_position[0], 2) + Math.pow(updated_position[1], 2) + Math.pow(updated_position[2], 2)) < radius_of_nucleus){
					sharp_positions[sharpMolecule] = updated_position;
				}			
			}
			final_num_sharp_inside = num_sharp_inside;
		}
				
		return (double)final_num_sharp_inside / (double)sharp_copy_number;
	}

	//run simulations for different ratios of diffusion coefficients
	private void runSimulations(double[] ratios){
		//ratios = [1,2,4,8,16,32,64,128,256,512,1024]
		for(int i=0; i<ratios.length; i++){
			runSimulation(ratios[i]);
		}
	}
	
}
