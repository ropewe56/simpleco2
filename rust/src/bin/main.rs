//use std::io;https://www.forrestthewoods.com/blog/how-to-debug-rust-with-visual-studio-code/

#![feature(non_ascii_idents)]
#![allow(non_snake_case)]
#![allow(unused)]

extern crate clap;
use clap::{Arg, App, SubCommand};
extern crate dirs;
use std::path::Path;
use std::path::PathBuf;
use std::env;

extern crate radtrans;
use radtrans::spectrum::spectrum;

fn main() { //  -> hdf5::Result<()>)

    // current working dir
    let cwd0 = env::current_dir().unwrap();
    let cwd  = cwd0.to_str().unwrap();

    // home dir
    let home_dir = dirs::home_dir().unwrap();
    let home_dir = home_dir.to_str().unwrap();

    //let jsonf : PathBuf = ["/home/wester/Projects/radtrans/radinput/", "input.json"].iter().collect();
    let jsonf : PathBuf = ["/home/rolf/Projects/Klima/radtrans/radoutput/b/0_1_1.0_1e-11/", "input.json"].iter().collect();
    let jsonf : PathBuf = ["/home/rolf/Projects/Klima/radtrans/radoutput/c/00_0_1_1.0_1e-11/", "input.json"].iter().collect();
    let mut json_file = jsonf.to_str().unwrap().to_string();

    for (i, argument) in env::args().enumerate() {
        println!("{}", argument);
        let vec: Vec<&str> = argument.split("=").collect();
        if vec[0] == "--json_file" {
            json_file = vec[1].to_string();
        }
    }

    let mut spec = radtrans::spectrum::spectrum::Spectrum::new(&json_file);
    spec.integrate();
    if spec.par.compute_F == true {
        spec.compute_F();
    }
}