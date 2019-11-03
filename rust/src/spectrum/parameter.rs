use serde_json::{Value};
use ndarray::prelude::*;
use ndarray::{Array};
use std::iter::FromIterator;

use std::io;
use std::io::prelude::*;
use std::fs::File;

pub fn read_json(jsonfname: &str) -> Result<Value, io::Error> {
    let mut f = File::open(jsonfname)?;
    let mut buffer = String::new();
    f.read_to_string(&mut buffer)?;
    let v: Value = serde_json::from_str(&buffer)?;
    Ok(v)
}

pub struct Parameter {
    pub h                   : f64,
    pub c0                  : f64,
    pub kB                  : f64,
    pub σ                   : f64,
    pub pi                  : f64,
    pub mCO2                : f64,
    pub lmin                : f64,
    pub lmax                : f64,
    pub dl                  : f64,
    pub surface_temperature : f64,
    pub temperature_mode    : i64,
    pub Dl_factor           : f64,
    pub p_ref               : f64,
    pub T_ref               : f64,
    pub albedo              : f64,
    pub with_emission       : i64,
    pub background          : f64,
    pub max_isotope_id      : usize,
    pub integrate           : bool,
    pub compute_F           : bool,
    
    pub out_dir        : String,
    pub logfile        : String,
    pub logfile_F      : String,
    pub infofile       : String,
    pub NCO2_path      : String,
    pub theta_path     : String,
    pub spdat_path     : String,
    pub T_Q_path       : String,
    pub h_p_path       : String,
    pub h_T_path       : String,
    pub z_path         : String,
    pub z_iout         : String,
}

fn get_bool(json: &Value, key1: &str, key2: &str) -> bool {
    let vv = &json[key1][key2];
    let vvf : Option<bool> = vv.as_bool();
    match vvf {
        None    => { panic!("key1 = {}, key2 = {}", key1, key2) },
        Some(v) => { v }
    }
}

fn get_i64(json: &Value, key1: &str, key2: &str) -> i64 {
    let vv = &json[key1][key2];
    let vvf : Option<i64> = vv.as_i64();
    match vvf {
        None    => { panic!("key1 = {}, key2 = {}", key1, key2) },
        Some(v) => { v }
    }
}

fn get_f64(json: &Value, key1: &str, key2: &str) -> f64 {
    let vv = &json[key1][key2];
    let vvf : Option<f64> = vv.as_f64();
    match vvf {
        None    => { panic!("key1 = {}, key2 = {}", key1, key2) },
        Some(v) => { v }
    }
}

fn get_f64_array(json: &Value, key1: &str, key2: &str) -> Array1<f64> {
    let vv = &json[key1][key2];
    let vvf : Option<&Vec<serde_json::value::Value>> = vv.as_array();
    match vvf {
        None    => { panic!("key1 = {}, key2 = {}", key1, key2) },
        Some(v) => { Array::from_iter( v.iter().map(|x| x.as_f64().unwrap()) ) }
    }
}

fn get_string(json: &Value, key1: &str, key2: &str) -> String {
    let vv = &json[key1][key2];
    let vvs :  Option<&str> = vv.as_str();
    match vvs {
        None    => { panic!("key1 = {}, key2 = {}", key1, key2) },
        Some(v) => { v.to_string() }
    }
}

impl Parameter {
    pub fn new(json_file: &str) -> Parameter {

        let mut par : Parameter;
        let json_ = read_json(&json_file);
        let json = match json_ {
            Ok(json) => { json }
            Err(err) => { println!("error: {}. Could not open {}", err, json_file); std::process::exit(1); }
        };

        Parameter{  h                   : get_f64(&json, "natconst", "h"),
                    c0                  : get_f64(&json, "natconst", "c0"),
                    kB                  : get_f64(&json, "natconst", "kB"),
                    σ                   : get_f64(&json, "natconst", "sigma"),
                    mCO2                : get_f64(&json, "natconst", "mCO2"),
                    pi                  : get_f64(&json, "natconst", "pi"),
                    lmin                : get_f64(&json, "run_parameter", "lmin"),
                    lmax                : get_f64(&json, "run_parameter", "lmax"),
                    dl                  : get_f64(&json, "run_parameter", "dl"),
                    Dl_factor           : get_f64(&json, "run_parameter", "Dl_factor"),
                    surface_temperature : get_f64(&json, "run_parameter", "surface_temperature"),
                    T_ref               : get_f64(&json, "run_parameter", "T_ref"),
                    p_ref               : get_f64(&json, "run_parameter", "p_ref"),
                    albedo              : get_f64(&json, "run_parameter", "albedo"),
                    temperature_mode    : get_i64(&json, "run_parameter", "temperature_mode"),
                    with_emission       : get_i64(&json, "run_parameter", "with_emission"),
                    background          : get_f64(&json, "run_parameter", "background"),
                    max_isotope_id      : get_i64(&json, "run_parameter", "max_isotope_id") as usize,
                    integrate           : get_bool(&json, "run_parameter", "integrate"),
                    compute_F           : get_bool(&json, "run_parameter", "compute_F"),
     
                    out_dir        : get_string(&json, "paths", "out_dir"),
                    logfile        : get_string(&json, "paths", "logfile"),
                    logfile_F      : get_string(&json, "paths", "logfile_F"),
                    infofile       : get_string(&json, "paths", "infofile"),
                    NCO2_path      : get_string(&json, "paths", "NCO2_path"),
                    theta_path     : get_string(&json, "paths", "theta_path"),
                    spdat_path     : get_string(&json, "paths", "spdat_path"),
                    T_Q_path       : get_string(&json, "paths", "T_Q_path"),
                    h_p_path       : get_string(&json, "paths", "h_p_path"),
                    h_T_path       : get_string(&json, "paths", "h_T_path"),
                    z_path         : get_string(&json, "paths", "z_path"),
                    z_iout         : get_string(&json, "paths", "z_iout"),
        }
    }
}
