#ifndef FILE_IO_H
#define FILE_IO_H

//------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>

//------------------------------------------------------------

namespace file_io_utils
{
  template <typename T> struct File
  {
    File() : v(false){};
    File(char* fname) : v(false) { open(fname); }
    File(bool verbose) : v(verbose) {} 
    File(char* fname, bool verbose) : v(verbose) { open(fname); }
    ~File() { if(fs.is_open()) close(); }
         
    T fs;
    bool v;
         
    void open(char* fname);
    void close(char* fname );
    void close();
  }; // struct File
   
  //------------------------------------------------------------
   
  template <typename T> void File<T>::open(char* fname)
  {
    if(v) std::cout << "... opening file: " << fname;
      
    try {
      fs.open(fname);
    }   
    catch (std::exception &e) {
      if(v)
        std::cerr << "... unable to open file " << fname << std::endl
                  << "... Standard exception: " << e.what() << std::endl;
      std::exit(1); 
    }
    if(v) std::cout << " ... done\n";
  }
   
  //------------------------------------------------------------
   
  template <typename T> void File<T>::close(char* fname)
  {
    if(v) std::cout << "... closing file: " << fname;
      
    try {
      fs.close();
    }      
    catch (std::exception &e) {
      if(v)      
        std::cerr << "... unable to close file " << fname << std::endl
                  << "... Standard exception: " << e.what() << std::endl;
      std::exit(1); 
    }   
    if(v) std::cout << " ... done\n";
  }      

  //------------------------------------------------------------
   
  template <typename T> void File<T>::close()
  {
    if(v) std::cout << "... closing file: " ;
      
    try {
      fs.close();
    }      
    catch (std::exception &e) {
      if(v)      
        std::cerr << "... unable to close file " << std::endl
                  << "... Standard exception: " << e.what() << std::endl;
      std::exit(1); 
    }
    if(v) std::cout << " ... done\n";
  }
} // namespace file_io_utils

//------------------------------------------------------------

#endif // FILE_IO_H


/* example 

namespace fio = file_io_utils;

int main(int argc, char* argv[])
{
std::string str;
char file_name[] = "tmp.dat";
bool verb =true;
   
{
      fio::File<std::ofstream> data;      
      data.open(file_name);
      data.fs << "akjsdhakljdfhakljdfh\n" ;   
      data.close(file_name);
   }
   
   {
      fio::File<std::ifstream> in;
      in.open(file_name);
      in.fs >> str;
      std::cout << str << std::endl;
      in.close(file_name);
   }

   {
      fio::File<std::ifstream> in2(file_name);
      //in.open(file_name);
      in2.fs >> str;
      std::cout << str << std::endl;
      in2.close(file_name);
   }

   {
      fio::File<std::ifstream> in3(file_name, verb);
      //in.open(file_name);
      in3.fs >> str;
      std::cout << str << std::endl;
      in3.close(file_name);
   }

   {
      fio::File<std::ifstream> in4(verb);
      in4.open(file_name);
      in4.fs >> str;
      std::cout << str << std::endl;
      in4.close(file_name);
   }

   return 0;
}

*/
