# SSTree

o/ 

A huachafo SSTree implementation for points. 

## Running

To run this project, need to adjust the CMakeFile the SFML directory.
```asm
set(SFML_DIR D:/SFML/lib/cmake/SFML)
```

You will also need to set the embedding.json route file on function int main()
```asm
const std::string FILE_NAME("../embedding.json");
```

Then compile and indexing will execute. Automatically it will create a file named embedding.dat (will take around 6 minutes). It is supossed to be there the SSTree. Adjust indexing.cpp as you wish :)


main.cpp was just for testing
