/* This game driver doing serial I/O with Game_project
 * Group 11 (Charlie Wang, Fredrick Song and Yizhou Luo) 
 * Mar 6, 2023
 * 
 * Example Sensor: Joystick
 *                         GND   connected to GND on Arduino
 *                         +5V   connected to +5V on Arduino
 *                         VRx   connected to  A0 on Arduino
 *                         VRy   connected to  A1 on Arduino
 *       Sensor: Button Switch                       
 *                               connected to Pin 7 on Arduino
 *                               
 */
const byte y_axis_PIN=A0;      // use A0 for joystick vertical-axis
const byte x_axis_PIN=A1;      // use A1 for joystick horizontal-axis
const byte switch_PIN1=7;     // use D7 for switch button
int i=1;                       // initialize Arduino counter variable to check for synchronization

void setup()  // The setup routine runs once when you press reset
{
  pinMode(x_axis_PIN,INPUT); // set joystick x-axis pin to be INPUT
  pinMode(y_axis_PIN,INPUT); // set joystick y-axis pin to be INPUT
  pinMode(switch_PIN1,INPUT); // set switch pin to be INPUT
  Serial.begin(115200);      // Initialize serial communication at 115200 bits per second
  delay(500);                // allow time for communication to be initialized
}

void loop() // The communication loop routine runs over and over again forever
{
  int x,y,s1; // x- and y-axis joystick values (0 to 1023, 511 middle position)
  int data_from_Matlab; // store counter from Matlab which can be used to check synchronization
  if(Serial.available()>2) // check if serial data has become available 
  {
    x=analogRead(x_axis_PIN); // read A0 x-axis joystick value
    y=analogRead(y_axis_PIN); // read A1 y-axis joystick value
    s1=digitalRead(switch_PIN1); // read D7 Switch button value as 0 or 1.
    // read incoming string data from Matlab consisting of counter to check for synchronization
    data_from_Matlab=Serial.parseInt();// load first valid integer in incoming serial
    Serial.print(String(String(i)+","+String(x)+","+String(y)+","+String(s1))); // send to Matlab Arduino counter to check for  // ","+String(s2)
                                                                 // synchronization, and joystick data
    Serial.write(13);          // "Carriage Return"
    Serial.write(10);          // "Linefeed"
    Serial.flush();            // wait until serial string is finished sending to Matlab
    i++;                       // increment Arduino counter to check for synchronization
  }
}
