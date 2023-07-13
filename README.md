#  ImgTool

A command line and/or scriptable image hackery tool.

RPN based, like old HP calculators.

Here's the current commands, obtained by running the tool without arguments.

```

Chris Lomont's RPN image tool v0.1
Usage: This is an RPN based image tool. Command args are RPN commands.
       Commands either on command line or run as --script filename
       --verbose to print more
       Each command shows what it does to the stack.
read        : filename -> image, loads image
write       : image filename -> ,  outputs saved image
image       : w h r g b a -> image, makes image size w x h, color rgba in 0-1
getpixel    : img i j -> img r g b a, reads pixel 0-1
setpixel    : img i j r g b a -> img, writes pixel 0-1
colorspace  : image space -> image', where space=[linear|sRGB|YCbCr|RGB], does conversion
error       : im1 im2 errtype -> im1 im2 errval, prints error type mse, psnr, ssim
maxc        : img -> max, max value of all r,g,b values in image
size        : img -> w h, where w,h is size in pixels
resize      : img w h style -> img', resize to w h by style nn,bilinear,bicubic,lanczos2,lanczos3,lanczos4
resize%     : img v style -> img', resize by v%, style as above
resize*     : img m style -> img', resize by multiplier m, style as above
gaussian    : img s -> img' , gaussian blur, std dev s
rotate      : TODO: img angle expand -> img', rotate image by angle degrees, expand true makes bigger to center, false keeps size
crop        : img x1 y1 x2 y2 -> img', crop image to rectangle (x1,y1)-(x2,y2) inclusive
pad         : img top bottom left right r g b a -> img2, pad image with given color, given pixel margins
flipx       : img -> img2, flip image
flipy       : img -> img2, flip image
files       : path regex -> f1 f2 ... fn n, reads files matching regex, pushes on stack with count
version     :  -> major minor, get version
timeus      :  -> time_us, get elapsed time in us
rand        : a b -> rand32(a,b), uniform random integer in [a,b)
srand       : seed -> , set random seed to integer seed
arg         :  n -> arg, get command line arg n, passed via -a item, n = 1,2,...
abs         : item -> abs(img)
ceil        : item -> ceil(item)
floor       : item -> floor(img)
round       : item -> round(item)
min         : a b -> min(a,b)
max         : a b -> max(a,b)
clamp       : item1 a b -> clamp(item1,a,b)
sin         : item -> abs(img)
cos         : item -> abs(img)
pi          :  -> pi
e           :   -> e
pow         : item1 item2 -> pow(item1,item2)
exp         : a -> e^a
log         : val base -> log_base(val)
neg         : a -> -a
sign        : a -> sign(a), is -1,0,1
+           : item1 item2 -> item1+item2
-           : item1 item2 -> item1-item2
*           : item1 item2 -> item1*item2
/           : item1 item2 -> item1/item2
mod         : a b -> a mod b
==          : item1 item2 -> item1==item2, 0 if false, else 1
!=          : item1 item2 -> item1!=item2, 0 if false, else 1
>=          : item1 item2 -> item1>=item2, 0 if false, else 1
<=          : item1 item2 -> item1<=item2, 0 if false, else 1
>           : item1 item2 -> item1>item2, 0 if false, else 1
<           : item1 item2 -> item1<item2, 0 if false, else 1
dup         : a -> a a, duplicates top item
dup2        : a b -> a b a b, duplicates top item
dupn        : x1 x2 .. xn n -> x1 x2 .. xn x1 x2 .. xn, duplicate top n
drop        : a -> , drops top item
drop2       : a b -> , drops top item
dropn       : x1 x2 .. xn n -> , drops top n
swap        : a b -> b a , swaps top 2 items
over        : a b -> a b a , copies item at level 2 to top
rot         : 3 2 1 -> 2 1 3 , rotates item in level 3 to level 1, 1 to 2, 2 to 3
unrot       : 3 2 1 -> 1 3 2 , rotates opposite of rot
roll        : x1 x2.. xn n -> x2 x3 ... xn x1  , like rot, but n items
rolld       : x1 x2 .. xn n -> xn x1 x2 x3 ... xn-1 , reverse of roll
pick        : xn ... x1 n -> xn ... x1 xn , copies item xn to top
depth       : ... -> n , pushes depth to top of stack
print       : item -> , prints top item
printn      : x1 x2 ... xn  n -> , prints top N items
endl        :  -> endline, pushes an endline string
label       : name -> , create named label for next item index
ja          : label -> , JumpAlways: goto label
je          : val label -> , JumpEqual: if val != 0, goto label
halt        :  exitcode -> , stops program, returns code
sto         : item name -> , store item in name
rcl         : name -> item, look up item
dumpstate   :  -> , print out state items
system      : cmd -> return_value, execute cmd on system call - WARNING - be careful!
verbosity   : v -> , set verbosity 0=none, 1=info, 2= all
if          : t1 t2 .. tn f1 f2 .. fm n m b -> ti or fj, if b != 0, keep t1..tn, else keep f1..fm
->str       : item -> 'item', formats item as string
rangeloop   : min max -> , loops over index in [min,max], each iter puts index on stack, use endloop
itemloop    : i1 i2 .. in n -> , loops over items in {i1,i2,..,in}, each iter puts item then index i=0+ on stack, use endloop
endloop     :  -> , ends loop, jumps to top
subroutine  :  name -> , starts subroutine, ends with endsub
endsub      :  name -> , ends subroutine
gosub       :  name -> , jumps to subroutine
return      :  -> , returns from subroutine

```