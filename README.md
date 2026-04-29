# MakeRotLib in Rust

A rotamer builder written in Rust for Rosetta which copies the MakeRotLib protocol written by P. Doug Renfrew and Andy Watkins.
This project is for my own edification on Rust and how Rosetta parameterizes & builds conformer libraries.
My plan is to make a another library after this one that improves upon the current methods available, but will be faster, more memory safe, and a standalone easy to use application.

Plan
---
- Work my way through rewriting all of the MakeRotLib protocol
- I will need to make an approximation of the Rosetta ScoreFunction as a way to score conformers

License
---
MIT - Feel free to use this for any project where you would like to make a library of param files to use in Rosetta.

---
Author: Andrew Powers (apowers4@uoregon.edu)
