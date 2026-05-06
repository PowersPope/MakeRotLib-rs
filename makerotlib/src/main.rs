/// @brief: A Rust recreation of the MakeRotLib protocol implemented by
/// Doug Renfrew and Andy Watkins originally
/// @author: Andrew Powers ( apowers4@uoregon.edu )

use tracing_subscriber;
use tracing::{info, Level};
pub mod makerotlib_options;

fn main() {
    // Init our tracing subscriber for trace info.
    let subscriber = tracing_subscriber::FmtSubscriber::builder()
        .with_max_level(Level::TRACE)
        .finish();
    tracing::subscriber::set_global_default(subscriber)
        .expect("Setting default subscriber failed");

    let mut out = makerotlib_options::read_in_data("inputs/C40_rot_lib_options_test.in");
    makerotlib_options::second_file_parse("inputs/C40_rot_lib_options_test.in", &mut out);
    info!("Our MakeRotLibOptionsData object is now constructed and ready to be used!");
}

