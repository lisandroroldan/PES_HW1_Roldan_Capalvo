#!/bin/sh -f


rm -f "$2/$1.boh"
rm -f "$2/$1.post.res"
rm -f "$2/$1.out"

# OutputFile: $2/$1.boh
# ErrorFile: $2/$1.err



"$3/frame3dd" "$2/$1.dat" "$2/$1.out" | tee "$2/$1.boh"




