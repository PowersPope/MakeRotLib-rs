/// @brief: A RotData struct type that will be used to load our data into when calling
/// the mover.
/// @author: Andrew Powers (apowers4@uoregon.edu)


pub struct RotData {
    num_bbs_:       u32,
    omega_ :        f32,
    min_omega_:     f32,
    epsilon_:       f32,
    energy_:        f32,
    probability_:   f32,
    num_chi_:       u32,
    num_clusters_:  u32,
    cluster_num_:   u32,
    twist_:         u32,
    inter_rep_:     f32,
    inter_atr_:     f32,
    intra_rep_:     f32,
    intra_atr_:     f32,
    solvation_:     f32,
    semirotameric_: bool,
    bbs_:           Vec<u32>,
    bb_ids_:        Vec<u32>,
    inp_chi_:       Vec<f32>, // Starting chi angles
    min_chi_:       Vec<f32>, // minimized chi angles
    lib_chi_val_:   Vec<u32>, // rotamer number for dunbrack format
    std_dev_:       Vec<f32>, // standard deviation of chi angles
    cen_dst_:       Vec<f32>, // distance from each centroid
}

