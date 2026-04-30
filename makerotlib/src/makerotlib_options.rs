/// @brief: Read in the options from a precomputed file to generate our parameters
/// that will be used when MakeRotLib is called
/// @author: Andrew Powers (apowers4@uoregon.edu)


use std::fs::OpenOptions;
use std::io::BufReader;
use std::io::prelude::*;

//// STRUCTS AND ENUM SETUP 

#[derive(Debug)]
/// Options that will be used to run
/// MakeRotLib, though this struct
/// needs to be loaded up with parameters
/// from a passed in file first
pub struct MakeRotLibOptions {
    n_bb_: u32,
    n_chi_: u32,
    n_centroids_: u32,
    semirotameric_: bool,
    kbt_: f64,
    name_: String,
    chi_ranges_: Vec<i32>,
    bb_ids_: Vec<i32>,
    bb_ranges_: Vec<i32>,
}

#[derive(Debug)]
// Describe range of torsion angle values
pub enum TorsionRange {
    Low(i32),
    High(i32),
    Step(i32),
}

pub enum CentroidRotNum {
    Angle(i32),
    RotNum(usize),
}

pub enum MakeRotLibPolymerType {
   PEPTIDE,
   PEPTOID,
}

////-------- FUNCTIONS BELOW ---------////

/// @brief: Take in a passed `@{filepath}` and convert the contents into
/// MakeRotLib parameters
/// @author: Andrew Powers
pub fn read_in_data( filepath: &str ) -> MakeRotLibOptions {

    let mut mklo = MakeRotLibOptions {
        n_bb_ : 0,
        n_chi_ : 0,
        n_centroids_ : 0,
        semirotameric_ : false,
        kbt_ : 0.0,
        name_ : "UNK".to_string(),
        chi_ranges_: Vec::new(),
        bb_ids_: Vec::new(),
        bb_ranges_: Vec::new(),
    };

    // Read in a passed file and error if it is not available
    let file_result =  OpenOptions::new().read(true).open(filepath);
    let infile = match file_result {
        Ok(infile) => infile,
        Err(e) => {
            panic!("failed to open file: {}", e.to_string());
        }
    };

    // Generate an iterable buffer object to iter through
    let reader = BufReader::new(infile);
    for l in reader.lines() {
        let line = l.expect("This is not a line that can be read.");

        if !line.starts_with("#") {
            let mut line_iter = line.split_whitespace();
            let tag: &str = line_iter.next()
                .expect("There is no next item or no whitespace to split");
            if tag == "NUM_CHI" {
                mklo.n_chi_ = line_iter.next()
                    .expect("No value after chi")
                    .parse::<u32>()
                    .unwrap();
            } else if tag == "NUM_BB" {
                mklo.n_bb_ = line_iter.next()
                    .expect("No value for BB")
                    .parse::<u32>()
                    .unwrap();
            } else if tag == "AA_NAME" {
                mklo.name_ = line_iter.next()
                    .expect("There is no name or no whitespace to split")
                    .to_string();
            } else if tag == "CENTROID" {
                mklo.n_centroids_ += 1;
            } else if tag == "TEMPERATURE" {
                mklo.kbt_ = line_iter.next()
                    .expect("No value for Temp")
                    .parse::<f64>()
                    .unwrap();
            } else if tag == "SEMIROTAMERIC" {
                mklo.semirotameric_ = true;
            } else if tag == "ROTAMERIC" {
                mklo.semirotameric_ = false;
            }
        }
    }

    // Update the sizes of our vectors, so that we can fill it with information
    mklo.chi_ranges_.resize(usize::try_from(mklo.n_chi_).unwrap(), 0);
    mklo.bb_ids_.resize(usize::try_from(mklo.n_bb_).unwrap(), 0);
    mklo.bb_ranges_.resize(usize::try_from(mklo.n_bb_).unwrap(), 0);

    println!("{:?}", mklo);
    return mklo
}
