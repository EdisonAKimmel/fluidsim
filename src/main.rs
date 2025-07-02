use std::{mem, usize, vec};
use std::io::{self, Write};
extern crate minifb;
use minifb::{Key, ScaleMode, Window, WindowOptions};
const SIZE: usize =100; //size of the simulation (with boundry cells)
const CELL_SIZE: usize = 10;
const SCREEN_SIZE: usize = (SIZE-2)*CELL_SIZE;
const DT: f32 = 0.05; //time step

fn add_source(dens_field: &mut Vec<f32>, source_field: &Vec<f32>) {
        dens_field.iter_mut()
        .zip(source_field.iter()) //.zip pairs each element of dens and source
        .for_each(|(dens, src)| *dens += src * DT); //adds to the density based on the src

    }
fn index(x: usize, y: usize)->usize{
        return (y)*(SIZE)+x;
    }
fn bindex(x: usize, y: usize)->usize{
        return (y)*(SCREEN_SIZE)+x;
    }
fn diffuse(boundry_type: i32, source_field: &mut Vec<f32>, prev_source_field: &Vec<f32>, diffusion_coefficient: f32) {
    let discretized_coefficient: f32 = DT * diffusion_coefficient * (((SIZE-2)*(SIZE-2)) as f32); //note that size is decreased by two to account for boundry   
    for _ in 0..19 { //we use 20 seidel solver steps for now
        for x in 1..SIZE-1 {
            for y in 1..SIZE-1 {
                let i = index(x,y);
                let i_up =  index(x,y+1);
                let i_down = index(x,y-1);
                let i_left = index(x-1,y);
                let i_right = index(x+1,y);
                let neighbors = source_field[i_up]+ source_field[i_down]+source_field[i_left]+source_field[i_right];
                source_field[i] = (prev_source_field[i]+discretized_coefficient*neighbors)/(1.0+4.0*discretized_coefficient);
            }
        }
        set_boundry(boundry_type,source_field);
    }
}
fn advect(prev_dens: &Vec<f32>, dens: &mut Vec<f32>, v_field: &(Vec<f32>,Vec<f32>)) {
    let dt_disc = DT * (SIZE as f32); // scale DT to grid units
    for x in 1..(SIZE - 1) {
        for y in 1..(SIZE - 1) {
            let i = index(x, y);
            // Backtrace position in the previous frame according to velocity field
            let x_prev = (x as f32) - v_field.0[i] * dt_disc;
            let y_prev = (y as f32) - v_field.1[i] * dt_disc;
            // Clamp to stay inside valid range [0.5, SIZE - 1.5]
            // Prevent sampling out of bounds during interpolation
            let x_prev = x_prev.clamp(0.5, (SIZE as f32) - 1.5);
            let y_prev = y_prev.clamp(0.5, (SIZE as f32) - 1.5);
            // Indices of the grid cell containing the backtraced position
            let x0 = x_prev.floor() as usize;
            let y0 = y_prev.floor() as usize;
            let x1 = x0 + 1;
            let y1 = y0 + 1;
            // Fractional part for interpolation
            let x_frac = x_prev - (x0 as f32);
            let y_frac = y_prev - (y0 as f32);
            // Fetch densities from four surrounding points in previous density field
            let bottom_left = prev_dens[index(x0, y0)];
            let bottom_right = prev_dens[index(x1, y0)];
            let top_left = prev_dens[index(x0, y1)];
            let top_right = prev_dens[index(x1, y1)];
            // Bilinear interpolation of density at backtraced position
            dens[i] =
            bottom_left * (1.0 - x_frac) * (1.0 - y_frac) +
            bottom_right * x_frac * (1.0 - y_frac) +
            top_left * (1.0 - x_frac) * y_frac +
            top_right * x_frac * y_frac;
        }
    }
}
fn project(v_field: &mut (Vec<f32>,Vec<f32>), pressure: &mut Vec<f32>, divergence: &mut  Vec<f32>) {
    let h= 1.0 / ((SIZE-2) as f32);
    for x in 1..SIZE-1 {
            for y in 1..SIZE-1 {
                let i = index(x,y);
                let i_up =  index(x,y+1);
                let i_down = index(x,y-1);
                let i_left = index(x-1,y);
                let i_right = index(x+1,y);
                divergence[i] = -0.5 * h * (v_field.0[i_right] - v_field.0[i_left] + v_field.1[i_up] - v_field.1[i_down]); //finds the divergence
                //using finite difference method
                pressure[i] = 0.0; //initializes pressure to zero
            }
        }
        set_boundry(0,divergence);
        set_boundry(0,pressure);
        for _ in 0..19 { //we use 20 seidel solver steps for now
        for x in 1..SIZE-1 {
            for y in 1..SIZE-1 {
                let i = index(x,y);
                let i_up =  index(x,y+1);
                let i_down = index(x,y-1);
                let i_left = index(x-1,y);
                let i_right = index(x+1,y);
                let neighbors = pressure[i_up]+ pressure[i_down]+pressure[i_left]+pressure[i_right];
                pressure[i] = (neighbors+divergence[i])/4.0;
            }
        }
        set_boundry(0,pressure);
    }
    for x in 1..SIZE-1 {
            for y in 1..SIZE-1 {
                let i = index(x,y);
                let i_up =  index(x,y+1);
                let i_down = index(x,y-1);
                let i_left = index(x-1,y);
                let i_right = index(x+1,y);
                v_field.0[i] -= 0.5 * (pressure[i_right]-pressure[i_left]) / h; //subtracts the gradient of pressure from the velocity
                v_field.1[i] -= 0.5 * (pressure[i_up]-pressure[i_down]) / h; //making it divergence free
            }
        }
    set_boundry(1, &mut v_field.0); //set the boundry conditions for the velocity 
    set_boundry(2, &mut v_field.1); //x and y components
}
fn set_boundry(boundry_type: i32, source_field: &mut Vec<f32>) {
    if boundry_type == 0 {
        for i in 1..(SIZE - 1) {
        // top and bottom edges
        source_field[0 * SIZE + i] = source_field[1 * SIZE + i]; // top edge
        source_field[(SIZE - 1) * SIZE + i] = source_field[(SIZE - 2) * SIZE + i]; // bottom edge
        // left and right edges
        source_field[i * SIZE + 0] = source_field[i * SIZE + 1]; // left edge
        source_field[i * SIZE + (SIZE - 1)] = source_field[i * SIZE + (SIZE - 2)]; // right edge
        //source_field[i * SIZE + (SIZE - 1)] = 0.0; //put a hole in the right side of the box
    }
    }
    if boundry_type == 1 { //for vector quantities (x)
        for i in 1..SIZE-1 {
            source_field[SIZE*i] = -source_field[SIZE*i+1]; //for left edge
            //source_field[SIZE*i] = 0.05; //put a wind on the left edge
            source_field[i * SIZE + (SIZE-1)] = -source_field[i * SIZE + (SIZE-2)];//for right edge
        }
    }
    if boundry_type == 2 { //for vector quantities (y)
        for i in 1..SIZE-1 {
            source_field[i] = -source_field[i+SIZE]; // for top edge
            source_field[(SIZE-1)*SIZE+i] = -source_field[(SIZE -2) * SIZE + i]; //for bottom edge
        }
    }
    //set the corners equal to average edge values
    source_field[0] = 0.5 * source_field[1] + source_field[SIZE + 1];
    source_field[SIZE - 1] = 0.5 * source_field[SIZE - 2] + source_field[SIZE * 2 - 2];
    source_field[SIZE * (SIZE - 1)] = 0.5 * source_field[SIZE * (SIZE - 1) + 1] + source_field[SIZE * (SIZE - 2)];
    source_field[SIZE * SIZE - 1] = 0.5 * source_field[SIZE * SIZE - 2] + source_field[SIZE * (SIZE - 2) + (SIZE - 1)];
    }
fn dec_vel(v_field: &mut (Vec<f32>, Vec<f32>)) {
    for vx in v_field.0.iter_mut() {
        *vx *= 0.99;
    }
    for vy in v_field.1.iter_mut() {
        *vy *= 0.99;
    }

}
fn dec_dens(dens: &mut Vec<f32>) {
    for i in dens.iter_mut() {
        *i *= 0.999;
    }
}
fn dens_step(prev_dens: &mut Vec<f32>, dens: &mut Vec<f32>, v_field: &(Vec<f32>,Vec<f32>), diffusion_constant: f32) {
    add_source(dens, &prev_dens);
    mem::swap(prev_dens,dens);
    diffuse(0, dens, &prev_dens, diffusion_constant);
    mem::swap(prev_dens,dens);
    advect(prev_dens, dens, v_field);
}
fn vel_step(prev_dens: &mut Vec<f32>, dens: &mut Vec<f32>, v_field: &mut (Vec<f32>, Vec<f32>), v_field_prev: &mut (Vec<f32>, Vec<f32>),  pressure: &mut Vec<f32>, divergence: &mut Vec<f32>, viscosity: f32) {
    add_source(&mut v_field.0, & v_field_prev.0);
    add_source(&mut v_field.1, & v_field_prev.1);
    mem::swap(v_field_prev, v_field);
    diffuse(1, & mut v_field.0, & v_field_prev.0, viscosity);
    diffuse(2, & mut v_field.1, & v_field_prev.1, viscosity);
    mem::swap(v_field_prev, v_field);
    project(v_field, pressure, divergence);
    mem::swap(v_field_prev, v_field);
    advect(prev_dens, dens, v_field);
    project(v_field, pressure, divergence);
    dec_vel(v_field); //decrease density and vel fields 
    dec_dens(dens);
    prev_dens.fill(0.0);
    v_field_prev.0.fill(0.0);
    v_field_prev.1.fill(0.0);
}

fn get_mouse_velocity(mouse_pos_prev: &(f32,f32), mouse_pos: &(f32,f32))->(f32,f32) {
    ((mouse_pos.0 - mouse_pos_prev.0)/DT,
    (mouse_pos.1 - mouse_pos_prev.1)/DT)
}
fn sim_to_screen2(input: &Vec<f32>) -> Vec<u32> {
    let mut output = vec!(0; (SCREEN_SIZE)*(SCREEN_SIZE));
    for y in 0..SCREEN_SIZE {
        for x in 0..SCREEN_SIZE {
            let padded_width = SIZE; // full size including padding
            let effective_cells = padded_width - 2; // exclude padding edges

            let sx = x as f32 / CELL_SIZE as f32;
            let sy = y as f32 / CELL_SIZE as f32;

            // Clamp so that ceil never accesses out-of-bounds (max floor is effective_cells - 1)
            let sx = sx.min((effective_cells - 1) as f32);
            let sy = sy.min((effective_cells - 1) as f32);

            let x0 = sx.floor() as usize + 1; // shift for padding
            let x1 = sx.ceil() as usize + 1;
            let y0 = sy.floor() as usize + 1;
            let y1 = sy.ceil() as usize + 1;

            let fx = sx - sx.floor();
            let fy = sy - sy.floor();
            
            let f00 = input[index(x0, y0)];
            let f10 = input[index(x1, y0)];
            let f01 = input[index(x0, y1)];
            let f11 = input[index(x1, y1)];

            let top = f00 * (1.0 - fx) + f10 * fx;
            let bottom = f01 * (1.0 - fx) + f11 * fx;

            let value = top * (1.0 - fy) + bottom * fy;
            output[bindex(x, y)] = value.round() as u32;

        }
    }       
    return output;
}
fn sim_to_screen(input: &Vec<f32>) -> Vec<u32> {
    let mut output = vec!(0; (SCREEN_SIZE)*(SCREEN_SIZE));
    for y in 0..SCREEN_SIZE {
        for x in 0..SCREEN_SIZE {
            let sx = ((x+CELL_SIZE) as f32 / CELL_SIZE as f32) as usize;
            let sy = ((y+CELL_SIZE) as f32 / CELL_SIZE as f32) as usize;
            let value = input[index(sx, sy)];
            output[bindex(x, y)] = value as u32;

        }
    }       
    return output;
}





fn main() {
    let mut mouse_pos_prev= (0.0f32,0.0f32);
    let mut mouse_pos= (0.0f32,0.0f32);
    let mut buffer = vec![0u32; (SCREEN_SIZE)*(SCREEN_SIZE)];


    let mut v_field = (vec![0.0; SIZE*SIZE],vec![0.0; SIZE*SIZE]);
    let mut v_field_prev = (vec![0.0; SIZE*SIZE],vec![0.0; SIZE*SIZE]);
    let mut dens = vec![100.0f32; SIZE*SIZE];
    let mut dens_prev = vec![0.0f32; SIZE*SIZE];
    let mut pressure = vec![0.0f32; SIZE*SIZE];
    let mut divergence = vec![0.0f32; SIZE*SIZE];

    let mut blerp = false;

    print!("Enable Bi-Linear Interpolation?\n");
    io::stdout().flush().unwrap(); // Make sure the prompt prints before input

    let mut input = String::new();
    io::stdin().read_line(&mut input).expect("Failed to read line");
    let input = input.trim(); // remove newline
    if input.to_ascii_lowercase() == "y" || input.to_ascii_lowercase() == "yes" {
        blerp = true;
    }


    let mut window = Window::new(
        "Fluid Sim, Made by Edison",
        SCREEN_SIZE,
        SCREEN_SIZE,
        WindowOptions {
            resize: false,
            scale_mode: ScaleMode::UpperLeft,
            ..WindowOptions::default()
        },
    )
    .expect("Unable to create the window");

    window.set_target_fps(400);

    let mut size = (0, 0);

    while window.is_open() && !window.is_key_down(Key::Escape) {
        if let Some((x,y))= window.get_mouse_pos(minifb::MouseMode::Discard) {
            std::mem::swap(&mut mouse_pos, &mut mouse_pos_prev);
            mouse_pos.0 = x/(CELL_SIZE as f32);
            mouse_pos.1 = y/(CELL_SIZE as f32);
        }
        if window.get_mouse_down(minifb::MouseButton::Left) {
            if let Some((x,y))= window.get_mouse_pos(minifb::MouseMode::Discard) { //adds the mouse velocity to the simulation velocity
                let sim_x = x/(CELL_SIZE as f32); //adjusts for the screen size
                let sim_y = y/(CELL_SIZE as f32);
                let click_vel_mult = 10.0 / (CELL_SIZE as f32);
                
                v_field_prev.0[index((sim_x) as usize, (sim_y) as usize)] = get_mouse_velocity(&mouse_pos_prev, &mouse_pos).0*click_vel_mult;
                v_field_prev.0[index((sim_x+1.0) as usize, (sim_y) as usize)] = get_mouse_velocity(&mouse_pos_prev, &mouse_pos).0*click_vel_mult;
                v_field_prev.0[index((sim_x-1.0) as usize, (sim_y) as usize)] = get_mouse_velocity(&mouse_pos_prev, &mouse_pos).0*click_vel_mult;
                v_field_prev.0[index((sim_x) as usize, (sim_y+1.0) as usize)] = get_mouse_velocity(&mouse_pos_prev, &mouse_pos).0*click_vel_mult;
                v_field_prev.0[index((sim_x) as usize, (sim_y-1.0) as usize)] = get_mouse_velocity(&mouse_pos_prev, &mouse_pos).0*click_vel_mult;

                v_field_prev.1[index((sim_x) as usize, (sim_y) as usize)] = get_mouse_velocity(&mouse_pos_prev, &mouse_pos).1*click_vel_mult;
                v_field_prev.1[index((sim_x+1.0) as usize, (sim_y) as usize)] = get_mouse_velocity(&mouse_pos_prev, &mouse_pos).1*click_vel_mult;
                v_field_prev.1[index((sim_x-1.0) as usize, (sim_y) as usize)] = get_mouse_velocity(&mouse_pos_prev, &mouse_pos).1*click_vel_mult;
                v_field_prev.1[index((sim_x) as usize, (sim_y+1.0) as usize)] = get_mouse_velocity(&mouse_pos_prev, &mouse_pos).1*click_vel_mult;
                v_field_prev.1[index((sim_x) as usize, (sim_y-1.0) as usize)] = get_mouse_velocity(&mouse_pos_prev, &mouse_pos).1*click_vel_mult;
                
                /* 
                v_field_prev.0[index((sim_x+1.0) as usize, (sim_y) as usize)] = 1.0*click_vel_mult;
                v_field_prev.0[index((sim_x-1.0) as usize, (sim_y) as usize)] = -1.0*click_vel_mult;

                v_field_prev.1[index((sim_x) as usize, (sim_y+1.0) as usize)] = 1.0*click_vel_mult;
                v_field_prev.1[index((sim_x) as usize, (sim_y-1.0) as usize)] = -1.0*click_vel_mult;
                */
                println!("The velocity is {},{} and the index is {}",v_field.0[index(sim_x as usize, sim_y as usize)],v_field.1[index(sim_x as usize, sim_y as usize)],index(sim_x as usize, sim_y as usize));
            }
        }
        if window.get_mouse_down(minifb::MouseButton::Right) {
           if let Some((x,y))= window.get_mouse_pos(minifb::MouseMode::Discard) { //adds fluid to the simulation
                let sim_x = (x/(CELL_SIZE as f32))+ 1.0 as f32; //adjusts for the screen size
                let sim_y = (y/(CELL_SIZE as f32)) + 1.0 as f32;
                let click_fluid = 1.5;
                if sim_y > 0 as f32 {dens_prev[index((sim_x) as usize, (sim_y-1.0) as usize)] = click_fluid/1.5;}
                if sim_y < SIZE as f32 {dens_prev[index((sim_x) as usize, (sim_y+1.0) as usize)] = click_fluid/1.5;}
                if sim_x < SIZE as f32 {dens_prev[index((sim_x+1.0) as usize, (sim_y) as usize)] = click_fluid/1.5; }
                if sim_x > 0 as f32 {dens_prev[index((sim_x-1.0) as usize, (sim_y) as usize)] = click_fluid/1.5;  }
                if sim_x > 0 as f32 && sim_y > 0 as f32 {dens_prev[index(sim_x as usize, sim_y as usize)] = click_fluid;}
                println!("The density is {} and the index is {}",dens[index(sim_x as usize, sim_y as usize)],index(sim_x as usize, sim_y as usize));
        }
        }
        let new_size = window.get_size();
        if new_size != size {
            size = new_size;
            buffer.resize(size.0 * size.1, 0);
        }
        if let Some((x,y))= window.get_mouse_pos(minifb::MouseMode::Discard) {
            let sim_x = x/(CELL_SIZE as f32); //adjusts for the screen size
            let sim_y = y/(CELL_SIZE as f32);
            //println!("The dens is {} and the index is {}",dens[index(sim_x as usize, sim_y as usize)],index(sim_x as usize, sim_y as usize));
        }
        let buffer: Vec<u32>;
        if blerp == true{
            buffer = sim_to_screen2(&(dens).iter().map(|&f| ((f.clamp(0.0, 1.0))*255.0)).collect::<Vec<f32>>());
        }
        else {
            buffer = sim_to_screen(&(dens).iter().map(|&f| ((f.clamp(0.0, 1.0))*255.0)).collect::<Vec<f32>>());
        }
        vel_step(&mut dens_prev, &mut dens, &mut v_field, &mut v_field_prev, &mut pressure, &mut divergence, 0.001);
        dens_step(&mut dens_prev, &mut dens, &v_field, 0.001);
        window
            .update_with_buffer(&buffer, new_size.0, new_size.1)
            .unwrap();
    }
}