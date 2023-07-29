# tristan-ff

Compilation:
```sh
# configure
cmake -B build -D user=<USERFILE>
#compile
cmake --build build -j
```

The executable will be located in `build/src/`.

> To compile in debug mode, add `-D DEBUG=ON` when configuring the code.
