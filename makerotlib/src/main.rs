/// @brief: A Rust recreation of the MakeRotLib protocol implemented by
/// Doug Renfrew and Andy Watkins originally
/// @author: Andrew Powers ( apowers4@uoregon.edu )

pub mod makerotlib_options;

fn main() {
    let _out: makerotlib_options::MakeRotLibOptions = makerotlib_options::read_in_data("inputs/C40_rot_lib_options_test.in");
}
