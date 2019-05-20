#\bin\bash
cd C;
if make; then
    gcc -g -Wall -o csteen src/*.c
    cd ..
else
    cd ..
    echo "Build failed"
    ret 1
fi
