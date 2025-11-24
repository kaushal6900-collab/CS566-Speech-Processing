


explaination:-

1.in this project i have used the function start record provided by the sir to record live audio
2.as previous project my logic is same as to read 2 words and store them in different arrays and then chech their zcr 
  to differentiate between yes and no 
3.but in this program there is a use of silent counter to differentiat gap between the 2 spoken words 
4. this silence counter will remain 0 until we are speaking something as long as we stop speaking this counter will 
   start increamenting.
5.by trial and error i have calculated how many frames i take pause between two words and set the threshould after 
  which the code will start recording second word 
6.as i am using the 'c' variable as to track which word is recording after recordig second word it will break the 
  loop
7. now our code will proccede tu calculate zcr of both words after calculation the if else condition will tell us
   that which word is yes