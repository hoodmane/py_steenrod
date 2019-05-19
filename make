#\bin\bash
cd C;
if make; then
    gcc -g -o test src/*.c
    cd ..
else
    cd ..
    echo "Build failed"
    ret 1
fi
