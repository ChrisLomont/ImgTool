#  ImgTool

A command line and/or scriptable image hackery tool.

Reverse Polish Notation (RPN) and stack based, like old HP calculators.

The input is a string of 'tokens'; those that are 'verbs' execute, those that are not are pushed onto a stack.

The stack can hold float64s, strings, and image pointers.

Each 'verb' pops 0 or more items on the stack, operates on them, then pushes 0 or more items back on the stack.

The list of tokens, after expanding all 'include's, is indexed by a program counter, allowing labels, loops, and other control flow.

Tokens can be given raw on the command line, or imported as scripts.

Good luck!

[Github](https://github.com/ChrisLomont/ImgTool)

## Image notes

No operations are performed on images without explicit direction. When images are loaded, they are in whatever the bytes in the image format meant, which is usually a non-linear, [sRGB](https://en.wikipedia.org/wiki/SRGB) color format. To do various operations correctly, you must explicitly do them. For many image processing tasks (filtering, resizing, rotation, [alpha compositing](https://en.wikipedia.org/wiki/Alpha_compositing)) you want to be a in a linear, pre-multiplied (also called associative) alpha format. To convert from sRGB, you must convert to a linear space first, **then** to a pre-multiplied alpha , then do your work. To write back out, you must reverse these in the correct manner: reverse the pre-multipled alpha, then convert to sRGB, then save. PNG is an excellent format for all of this.

Note that alpha in files is treated as linear, even if color components are not, which is standard. PNG requires that RGBA alpha is linear, and colors are **not** pre-multiplied alpha.

Currently the PNG loader is [stb_image](https://github.com/nothings/stb/blob/master/stb_image.h), which does not support other color spaces. This may be changed in the future. You can always write your own converters.

## Commands

Here's the current commands, obtained by running the tool without arguments.

```
Chris Lomont's RPN image tool v0.5, https://github.com/ChrisLomont/ImgTool
Usage: This is an RPN based image tool. Command args are RPN commands.
       Commands either on command line or run as -s filename
       --verbose to print more, 0=none, 1=info, 2=all
       Each command shows what it does to the stack.
read        : filename -> image, loads image
write       : image filename -> ,  outputs saved image
image       : w h r g b a -> image, makes image size w x h, color rgba in 0-1
getpixel    : img i j -> img r g b a, reads pixel 0-1
setpixel    : img i j r g b a -> img, writes pixel 0-1
colorspace  : image space -> image', where space=[linear|sRGB|YCbCr|RGB], does conversion
alpha*      : image -> image', applies alpha pre-multiplication
alpha/      : image -> image', reverses alpha pre-multiplication (0 alpha -> 0,0,0,0 color)
error       : im1 im2 errtype -> im1 im2 errval, prints error, errtype mse, psnr, ssim
maxc        : img -> max, max value of all r,g,b values in image
size        : img -> w h, where w,h is size in pixels
resize      : img w h style -> img', resize to w h by style nn,bilinear[+],bicubic[+],lanczos2,lanczos3[+],lanczos4,lanczos2r,lanczos3r,lanczos4r
resize%     : img v style -> img', resize by v%, style as above
resize*     : img m style -> img', resize by multiplier m, style as above
gaussian    : img radius -> img' , gaussian blur with given radius
rotate      : img angle filter -> img', rotate image by angle degrees using filter shear3, nn, bilinear, bicubic, dct
shift       : img dx dy filter -> img', shift image by dx dy using filter (todo all nn for now)
crop        : img x1 y1 x2 y2 -> img', crop image to rectangle (x1,y1)-(x2,y2) inclusive
pad         : img top bottom left right r g b a -> img2, pad image with given color, given pixel margins
flipx       : img -> img2, flip image
flipy       : img -> img2, flip image
blit        : src dst -> dst', copy pixels from src to dst
blitc       : src dst dx dy -> dst' copy src pixels to dst, placing dest corner at dx dy
blitr       : src x1 y1 w h dst dx dy -> dst', copy rect from src x1 y1 w h to dst at dx dy
blitover    : src dst -> dst dx dy', alpha blend src OVER dst, at dx dy
boundary    : img [r g b a] mode -> img', set sample boundary mode to color (with rgba), clamp, reflect, reverse, tile
f->i        : f1 f2 .. fn n -> i1 i2 .. in, converts n values in 0-1 to n values in 0-255, useful for colors
i->f        : i1 i2 .. in n -> f1 f2 .. fn, converts n values in 0-255 to n values in 0-1, useful for colors
line        : img x1 y1 x2 y2 r g b a -> img with line
circle      : img x1 y1 radius r g b a -> img with circle
circlef     : img x1 y1 radius r g b a -> img with filled circle
rect        : img x1 y1 x2 y2 r g b a -> img with rectangle
rectf       : img x1 y1 x2 y2 r g b a -> img with filled rectangle
text        : img x1 y1 r g b a text 0 m -> img x2 y2, draws text in font (always 0), pixel size m, returns img and final position
csvstart    :  csvname header1 header2 ... n -> , start a CSV file with given headers
csvput      :  val header csvname -> , stores val under header name in named csv
csvwrite    :  csvname filename -> ,
abs         : item -> abs(img)
ceil        : item -> ceil(item)
floor       : item -> floor(img)
round       : item -> round(item)
min         : a b -> min(a,b)
max         : a b -> max(a,b)
clamp       : item a b -> clamp(item,a,b)
sin         : item -> sin(item), values in radians
cos         : item -> cos(item), values in radians
pi          :  -> pi
e           :   -> e
pow         : item1 item2 -> pow(item1,item2)
sqrt        : item1 -> sqrt(item)
exp         : item -> e^item
log         : val base -> log_base(val)
neg         : item -> -item
sign        : item -> sign(item), is -1,0,1
+           : item1 item2 -> item1 + item2
-           : item1 item2 -> item1 - item2
*           : item1 item2 -> item1 * item2
/           : item1 item2 -> item1 / item2
mod         : item1 item2 -> item1 mod item2
==          : item1 item2 -> item1 == item2, 0 if false, else 1
!=          : item1 item2 -> item1 != item2, 0 if false, else 1
>=          : item1 item2 -> item1 >= item2, 0 if false, else 1
<=          : item1 item2 -> item1 <= item2, 0 if false, else 1
>           : item1 item2 -> item1 > item2, 0 if false, else 1
<           : item1 item2 -> item1 < item2, 0 if false, else 1
and         : a b -> (a and b), bitwise 'and' on integers
or          : a b -> (a or b), bitwise 'or' on integers
xor         : a b -> (a xor b), bitwise 'xor' on integers
not         : a -> (not a), treating 0 as false, != 0 as true, boolean not
files       : path regex -> f1 f2 ... fn n, reads files matching regex, pushes on stack with count
version     :  -> major minor, get version
timeus      :  -> time_us, get elapsed time in us
rand        : a b -> rand32(a,b), uniform random integer in [a,b)
srand       : seed -> , set random seed to integer seed
arg         :  n -> arg, get command line arg n, passed via -a item, n = 1,2,...
argcount    :   -> argcount, count of command line args passed via via -a
include     :  filename -> , include file as text, each file included at most once
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
endl        : -> endline, pushes an endline string
label       : name -> , create named label for next item index
ja          : label -> , JumpAlways: goto label
je          : item1 item2 label -> , jump to label if item1 == item2
jne         : item1 item2 label -> , jump to label if item1 != item2
jlt         : item1 item2 label -> , jump to label if item1 < item2
jle         : item1 item2 label -> , jump to label if item1 <= item2
jgt         : item1 item2 label -> , jump to label if item1 > item2
jge         : item1 item2 label -> , jump to label if item1 >= item2
halt        : exitcode -> , stops program, returns code
sto         : item name -> , store item in name
rcl         : name -> item, look up item
dumpstate   :  -> , print out state items
system      : cmd -> return_value, execute cmd on system call - WARNING - be careful!
verbosity   : v -> , set verbosity 0=none, 1=info, 2=all
if          : t1 t2 .. tn f1 f2 .. fm n m b -> ti or fj, if b != 0, keep t1..tn, else keep f1..fm
->str       : item -> 'item', formats item as string
rangeloop   : min max -> , loops over index in [min,max], each iter puts index on stack, use endloop
itemloop    : i1 i2 .. in n -> , loops over items in {i1,i2,..,in}, each iter puts item then index i=0+ on stack, use endloop
endloop     :  -> , ends loop, jumps to top
subroutine  :  name -> , starts subroutine, ends with endsub
endsub      :  name -> , ends subroutine
gosub       :  name -> , jumps to subroutine
return      :  -> , returns from subroutine
->list      :  item1 item2 ... itemn n  -> list of items, convert n items into a list
list->      :  list -> item1 item2 ... itemn n, list of n items out
listlen     :  list -> list list_length, get length of list
listget     :  list k -> list item_k, get kth item from list, 0 indexed
listset     :  list item k -> list, set kth item from list, 0 indexed
sublist     :  list a b -> sublist, get sublist of items a (inclusive) to b (exclusive) 0 indexed
listins     :  list item k -> list, insert item at index k, 0 indexed
listdel     :  list k -> list, delete the k item, 0 indexed
listappend  :  list item -> append item to list
listjoin    :  list1 list2 -> list, join lists 1 and 2
```