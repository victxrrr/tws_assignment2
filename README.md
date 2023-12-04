<h1 align="center">
  <br>
  Scientific Software / Technisch-wetenschappelijke software
  <br>
</h1>

<h4 align="center">Assignment 2</h4>

<p align="center">
  <a href="#usage">Usage</a> •
  <a href="#author">Author</a>
</p>

## Usage

The following commands are available.

```bash
# To run the program corresponding to the second question of this assignment 
$ make q2
$ ./q2.out [N] [T]

# To run the program corresponding to the third question of this assignment
$ make q3
$ ./q3.out [N] [T] [delta_beta] [target]

# To run the program corresponding to the fourth question of this assignment
$ make q4
$ ./q4.out [filename]

# To generate a plot in pdf format of the last simulation ran with ./q2.out 
$ make plot

# To clean the directory
$ make clean
```
> **Note**
> Do not forget to `make clean` after switching the floating point precision in `utils.f90`.


## Author

Fully implemented by Victor Lepère.

---

> [Victor Lepère](mailto:victor.lepere@student.uclouvain.be) &nbsp;&middot;&nbsp;
> GitHub [@victxrrr](https://github.com/victxrrr)