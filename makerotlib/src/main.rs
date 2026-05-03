/// @brief: A Rust recreation of the MakeRotLib protocol implemented by
/// Doug Renfrew and Andy Watkins originally
/// @author: Andrew Powers ( apowers4@uoregon.edu )

pub mod makerotlib_options;

fn main() {
    let (mut out, mut rotwells_nchi) = makerotlib_options::read_in_data("inputs/C40_rot_lib_options_test.in");
    makerotlib_options::second_file_parse("inputs/C40_rot_lib_options_test.in", &mut out, &mut rotwells_nchi);

    println!("---Inside of main---");
    println!("{:?}", out);
    println!("{rotwells_nchi:?}");
}

