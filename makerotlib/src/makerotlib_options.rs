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
pub struct MakeRotLibOptionsData {
    n_bb_: u32,
    n_chi_: u32,
    n_centroids_: u32,
    semirotameric_: bool,
    kbt_: f64,
    name_: String,
    chi_ranges_: Vec<TorsionRange>,
    bb_ids_: Vec<i32>,
    bb_ranges_: Vec<TorsionRange>,
    omg_range_: TorsionRange,
    eps_range_: TorsionRange,
    centroid_data_: Vec<Vec<CentroidRotNum>>,
    polymer_type_: MakeRotLibPolymerType,
}

#[derive(Clone,Debug)]
// Describe range of torsion angle values
pub struct TorsionRange {
    low: i32,
    high: i32,
    step: i32,
}

impl TorsionRange {
    fn new(low: i32, high: i32, step: i32) -> Self {
        Self { low, high, step }
    }
}

#[derive(Clone,Debug)]
pub struct CentroidRotNum {
    angle: i32,
    rot_num: usize,
}

impl CentroidRotNum {
    fn new(angle: i32, rot_num: usize) -> Self {
        Self { angle, rot_num }
    }
}


#[derive(Debug)]
pub enum MakeRotLibPolymerType {
    UNKNOWN,
    PEPTIDE,
    PEPTOID,
}

////-------- FUNCTIONS BELOW ---------////

/// @brief determine the number of rotamer chi angles based on
/// if the amino acid is `@{semirot}` the number of chi rotamers are
/// based on `@{n_chi}`.
pub fn determine_nrotchi( semirot: bool, n_chi: &u32 ) -> usize {
    let nrotchi: usize = if semirot {
        usize::try_from(*n_chi).expect("n_chi_ is not u32") - 1
    } else {
        usize::try_from(*n_chi).expect("n_chi_ is not u32")
    };
    nrotchi
}

/// @brief: Take in a passed `@{filepath}` and convert the contents into
/// MakeRotLib parameters
/// @author: Andrew Powers
pub fn read_in_data( filepath: &str ) -> MakeRotLibOptionsData {

    let mut mklo = MakeRotLibOptionsData {
        n_bb_ : 0,
        n_chi_ : 0,
        n_centroids_ : 0,
        semirotameric_ : false,
        kbt_ : 0.0,
        name_ : "UNK".to_string(),
        chi_ranges_: vec![TorsionRange::new(0,0,0)],
        bb_ids_: Vec::new(),
        bb_ranges_: vec![TorsionRange::new(0,0,0)],
        omg_range_: TorsionRange::new(0,0,0),
        eps_range_: TorsionRange::new(0,0,0),
        centroid_data_: vec![vec![CentroidRotNum::new(0,0)]],
        polymer_type_: MakeRotLibPolymerType::UNKNOWN,
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
    mklo.bb_ids_.resize(usize::try_from(mklo.n_bb_).unwrap(), 0);
    mklo.chi_ranges_.resize(usize::try_from(mklo.n_chi_).unwrap(), TorsionRange::new(0,0,0));
    mklo.bb_ranges_.resize(usize::try_from(mklo.n_bb_).unwrap(), TorsionRange::new(0,0,0));

    println!("{:?}", mklo);
    return mklo
}

/// @brief: Now that our `@{mklo}` dataset has been setup and sizes specified, we can now
/// loop back through `@{filepath}` to grab out the information that we skipped.
/// @author: Andrew Powers
pub fn second_file_parse( filepath: &str, mklo: &mut MakeRotLibOptionsData ) {
    let mut bb_i: usize = 0;
    let nrotchi: usize = determine_nrotchi( mklo.semirotameric_, &mklo.n_chi_);
    let mut rotwells_for_chi: Vec<Vec<i32>> = vec![vec![0]];
//     rotwells_for_chi.resize(nrotchi, vec![0]);
    let mut rotwells_specified: bool = false;

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
            if tag == "PHI_RANGE" {
                mklo.bb_ranges_[bb_i].low = line_iter.next()
                    .expect("PHI RANGE LOW")
                    .parse::<i32>()
                    .expect("Can't unpack PHI RANGE LOW to i32");
                mklo.bb_ranges_[bb_i].high = line_iter.next()
                    .expect("PHI RANGE HIGH")
                    .parse::<i32>()
                    .expect("Can't unpack PHI RANGE HIGH to i32");
                mklo.bb_ranges_[bb_i].step = line_iter.next()
                    .expect("PHI RANGE STEP")
                    .parse::<i32>()
                    .expect("Can't unpack PHI RANGE STEP to i32");
                mklo.bb_ids_[bb_i] = 2;
                bb_i += 1;
            } else if tag == "PSI_RANGE" {
                mklo.bb_ranges_[bb_i].low = line_iter.next()
                    .expect("PHI RANGE LOW")
                    .parse::<i32>()
                    .expect("Can't unpack PHI RANGE LOW to i32");
                mklo.bb_ranges_[bb_i].high = line_iter.next()
                    .expect("PHI RANGE HIGH")
                    .parse::<i32>()
                    .expect("Can't unpack PHI RANGE HIGH to i32");
                mklo.bb_ranges_[bb_i].step = line_iter.next()
                    .expect("PHI RANGE STEP")
                    .parse::<i32>()
                    .expect("Can't unpack PHI RANGE STEP to i32");
                mklo.bb_ids_[bb_i] = 3;
                bb_i += 1;
            } else if tag == "BB_RANGE" {
                mklo.bb_ranges_[bb_i].low = line_iter.next()
                    .expect("BB RANGE LOW")
                    .parse::<i32>()
                    .expect("Can't unpack BB RANGE LOW to i32");
                mklo.bb_ranges_[bb_i].high = line_iter.next()
                    .expect("BB RANGE HIGH")
                    .parse::<i32>()
                    .expect("Can't unpack BB RANGE HIGH to i32");
                mklo.bb_ranges_[bb_i].step = line_iter.next()
                    .expect("BB RANGE STEP")
                    .parse::<i32>()
                    .expect("Can't unpack BB RANGE STEP to i32");
                bb_i += 1;
            } else if tag == "OMG_RANGE" {
                mklo.omg_range_.low = line_iter.next()
                    .expect("BB RANGE LOW")
                    .parse::<i32>()
                    .expect("Can't unpack BB RANGE LOW to i32");
                mklo.omg_range_.high = line_iter.next()
                    .expect("BB RANGE HIGH")
                    .parse::<i32>()
                    .expect("Can't unpack BB RANGE HIGH to i32");
                mklo.omg_range_.step = line_iter.next()
                    .expect("BB RANGE STEP")
                    .parse::<i32>()
                    .expect("Can't unpack BB RANGE STEP to i32");
            } else if tag == "EPS_RANGE" {
                mklo.eps_range_.low = line_iter.next()
                    .expect("BB RANGE LOW")
                    .parse::<i32>()
                    .expect("Can't unpack BB RANGE LOW to i32");
                mklo.eps_range_.high = line_iter.next()
                    .expect("BB RANGE HIGH")
                    .parse::<i32>()
                    .expect("Can't unpack BB RANGE HIGH to i32");
                mklo.eps_range_.step = line_iter.next()
                    .expect("BB RANGE STEP")
                    .parse::<i32>()
                    .expect("Can't unpack BB RANGE STEP to i32");
            } else if tag == "CHI_RANGE" {
                // get chi index
                let chi_num: usize = line_iter.next()
                    .expect("There is no next value for chi_num")
                    .parse::<usize>()
                    .expect("chi_num cant be usize") - 1;
                mklo.chi_ranges_[chi_num].low = line_iter.next()
                    .expect("CHI RANGE LOW")
                    .parse::<i32>()
                    .expect("Can't unpack CHI RANGE LOW to i32");
                mklo.chi_ranges_[chi_num].high = line_iter.next()
                    .expect("CHI RANGE HIGH")
                    .parse::<i32>()
                    .expect("Can't unpack CHI RANGE HIGH to i32");
                mklo.chi_ranges_[chi_num].step = line_iter.next()
                    .expect("CHI RANGE STEP")
                    .parse::<i32>()
                    .expect("Can't unpack CHI RANGE STEP to i32");
            } else if tag == "CENTROID" {
                let mut temp_crnv: Vec<CentroidRotNum> = vec![CentroidRotNum::new(0, 0)];
                if mklo.semirotameric_ {
                   temp_crnv.resize( 
                       usize::try_from(mklo.n_chi_).expect("n_chi_ is not u32") - 1, 
                       CentroidRotNum::new(0, 0),
                       );
                   for i in 0..usize::try_from(mklo.n_chi_).unwrap()-1 {
                    temp_crnv[i].angle = line_iter.next().expect("No Centroid ANGLE")
                        .parse::<i32>().expect("Cant unpack Centroid Angle to i32");
                    temp_crnv[i].rot_num = line_iter.next()
                        .expect("No Centroid rot_num")
                        .parse::<usize>().expect("Cant unpack Centroid rot_num to usize");
                   }
                } else {
                   temp_crnv.resize( 
                       usize::try_from(mklo.n_chi_).expect("n_chi_ is not u32"), 
                       CentroidRotNum::new(0, 0),
                       );
                   for i in 0..usize::try_from(mklo.n_chi_).unwrap() {
                       temp_crnv[i].angle = line_iter.next().expect("No Centroid ANGLE")
                           .parse::<i32>().expect("Cant unpack Centroid Angle to i32");
                       temp_crnv[i].rot_num = line_iter.next()
                           .expect("No Centroid rot_num")
                           .parse::<usize>().expect("Cant unpack Centroid rot_num to usize");
                   }
                }
                mklo.centroid_data_.push(temp_crnv);
            } else if tag == "ROTWELLS" {
                rotwells_specified = true;
                rotwells_for_chi.resize(nrotchi, vec![0]);
                let chi_num: usize = line_iter.next().expect("No chi_num")
                    .parse::<usize>().expect("Cant convert chi_num to uszie");
                let n_rotwells: usize = line_iter.next().expect("No rotwells")
                    .parse::<usize>().expect("Cant convert rotwells to usize");
                for _ in 0..n_rotwells {
                    let rotwell_angle: i32 = line_iter.next()
                        .expect("rotwell isnt correct")
                        .parse::<i32>().expect("Cant convert rotwell to i32");
                    rotwells_for_chi[chi_num].push(rotwell_angle);
                }
            }
        }
    }

    // Now make sure we selected either centroid or rotamer wells, not both
    assert_ne!( mklo.centroid_data_.len() > 0, rotwells_for_chi.len() > 1,
        "Warning: you specified both centroids and rotamer well combinations. We can't do both.");
    let nrotchi: usize = determine_nrotchi( mklo.semirotameric_, &mklo.n_chi_);

    // Generate all combinations of rotamer wells.
    // indices is a vector over all chi that says which rotwell we are workign with
    // for each chi
    if rotwells_specified {
        let mut indices: Vec<usize> = Vec::with_capacity(nrotchi);
        indices.resize(nrotchi, 0);
        rotwells_for_chi.push( vec![1,0] );

        let mut p: usize = 0;
        while indices[nrotchi] == 0 {
            let mut temp_crnv: Vec<CentroidRotNum> = vec![CentroidRotNum::new(0, 0)];
            temp_crnv.resize(nrotchi, CentroidRotNum::new(0,0));
            for i in 0..nrotchi {
                temp_crnv[i].angle = rotwells_for_chi[i][indices[i]];
                temp_crnv[i].rot_num = indices[i];
            }
            mklo.centroid_data_.push(temp_crnv);
            indices[0] += 1;
            while indices[p] > rotwells_for_chi[p].len() {
                if p < nrotchi {indices[p] = 1;}
                else { indices[p] = 0;}
                p = p + 1;
                indices[p] += 1;
                if indices[p] < rotwells_for_chi[p].len() {p=1;}
            }
        }
    }

    // Currently leaving out the Chemical get residue name_map type
    // As I have not implemented the ChemicalManager into Rust or ResidueType
    // for now we will HARDCODE a polyer_type_ of PEPTIDE.
    mklo.polymer_type_ = MakeRotLibPolymerType::PEPTIDE;
}
