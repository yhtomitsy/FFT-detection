#include "nrf.h"
#include "math.h"

#ifdef NRF51

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include "nrf_adc.h"
#include "boards.h"
#include "app_uart.h"
#include "app_error.h"
#include "nrf_delay.h"

#define MAX_TEST_DATA_BYTES     (15U)                /**< max number of test bytes to be used for tx and rx. */
#define UART_TX_BUF_SIZE 256                         /**< UART TX buffer size. */
#define UART_RX_BUF_SIZE 1                           /**< UART RX buffer size. */

#ifndef NRF_APP_PRIORITY_HIGH
#define NRF_APP_PRIORITY_HIGH                   1
#endif

//// config FFT parameters
#define ARRAY_SIZE                              64                              // the size of FFT
#define LOGN_ARRAY_SIZE                         6
#define BIN_STEP                                50000/ARRAY_SIZE/2
//// ATTENTION! It's important to configure fft offset below. 
//// This is the adc value you get in total silence in float. 
//// My mic have few volts offset do this value is 172
#define FFT_OFFSET                              84.0                           // value to distract 
//// Cry detection sensivity parameters
#define AVER_POWER                              4                               // average power of signal shoulâ be higher
#define HAPR                                    25                              // harmonic to average power ratio of signal should be higher
#define LOW_BAND                                700                             // pitch freq of signal should be between 700 Hz ...
#define HIGH_BAND                               4000                            // ... and 4000 Hz
#define SILENCE_TRSHLD                          5                               // signal with power lower than this considered as silence
#define NUM_SILENT_SAMPLES                      15                              // if it's a baby's cry there should more than 15 silent samples

volatile bool data1Ready = false, data2Ready = false;
float inArrayR[2][ARRAY_SIZE] = {0};                                            // vars to process FFT 
float inArrayI[ARRAY_SIZE] = {0};                                               // vars to process FFT 
float offset = 172.0;

void calcNewOffset(float * offset) {
  float val = 0.0;
  for (uint8_t i = 0; i < ARRAY_SIZE; i++) 
    val += inArrayR[0][i] + *offset;
  val /= ARRAY_SIZE;
  *offset = (*offset + val) / 2;
}       

// ADC handler fills 2 arrays with data
void ADC_IRQHandler(void) {
  static uint32_t sampleCount = 0, measureArr = 1;
  //offset = 172.0;          // auto offset calculation
  nrf_adc_conversion_event_clean();
  
  if (data1Ready == false || data2Ready == false) {
    /*uint16_t adcResult = nrf_adc_result_get();
    if (adcResult < FFT_OFFSET) adcResult = 0;
    else adcResult = adcResult - FFT_OFFSET;
    inArrayR[0][sampleCount] = (float)adcResult;*/

    inArrayR[0][sampleCount] = (float) nrf_adc_result_get() - offset;
    sampleCount++;
    
    if (sampleCount >= ARRAY_SIZE && measureArr == 1) {
      data1Ready = true;                // first frame ready to analyze
      measureArr = 2;
      calcNewOffset(&offset); 
    }
    if (sampleCount >= ARRAY_SIZE * 2 && measureArr == 2) {
      data2Ready = true;                // second frame ready to analyze
      sampleCount = 0;
      measureArr = 1;
      calcNewOffset(&offset); 
    }
  }    
  nrf_adc_start();                      // trigger next ADC conversion
}

void adcConfig(void) {
    const nrf_adc_config_t nrf_adc_config = { 
      NRF_ADC_CONFIG_RES_8BIT,                    // 20 us ack time -> 50kHz sample rate
      NRF_ADC_CONFIG_SCALING_INPUT_TWO_THIRDS, 
      NRF_ADC_CONFIG_REF_SUPPLY_ONE_HALF 
    };
    nrf_adc_configure( (nrf_adc_config_t *)&nrf_adc_config);
    nrf_adc_input_select(NRF_ADC_CONFIG_INPUT_7);
    nrf_adc_int_enable(ADC_INTENSET_END_Enabled << ADC_INTENSET_END_Pos);
    NVIC_SetPriority(ADC_IRQn, NRF_APP_PRIORITY_HIGH);
    NVIC_EnableIRQ(ADC_IRQn);
}


// NAME:          FFT.
// PARAMETERS:  
//    float *Rdat    [in, out] - Real part of Input and Output Data (Signal or Spectrum)
//    float *Idat    [in, out] - Imaginary part of Input and Output Data (Signal or Spectrum)
//    int    N       [in]      - Input and Output Data length (Number of samples in arrays)
//    int    LogN    [in]      - Logarithm2(N)
// RETURN VALUE:  false on parameter error, true on success.
//_________________________________________________________________________________________
// NOTE: In this algorithm N and LogN can be only:
//       N    = 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384;
//       LogN = 2, 3,  4,  5,  6,   7,   8,   9,   10,   11,   12,   13,    14;
//_________________________________________________________________________________________
// FFT & ADC CONFIGURATION:
// sample rate - 50000 Hz (from ADC configuration)
// bins - 50000/ARRAY_SIZE/2 - 195,3125 Hz if ARRAY_SIZE = 128
// bin1: -97,65625..97,65625  Hz (center 0 Hz)
// bin2: 97,65625 ..292,96875 Hz (center 195,3125)
// bin3: 292,96875..488,28125 Hz
// bin4: 488,28125..683,59375 Hz
// bin5...
#define  NUMBER_IS_2_POW_K(x)   ((!((x)&((x)-1)))&&((x)>1))  // x is pow(2, k), k=1,2, ...

bool  FFT(float *Rdat, float *Idat, int N, int LogN) {
  // parameters error check:
  /*if((Rdat == NULL) || (Idat == NULL))                  return false;
  if((N > 16384) || (N < 1))                            return false;
  if(!NUMBER_IS_2_POW_K(N))                             return false;
  if((LogN < 2) || (LogN > 14))                         return false;*/

  register int  i, j, n, k, io, ie, in, nn;
  float         ru, iu, rtp, itp, rtq, itq, rw, iw, sr;
  
  static const float Rcoef[14] =
  {  -1.0000000000000000F,  0.0000000000000000F,  0.7071067811865475F,
      0.9238795325112867F,  0.9807852804032304F,  0.9951847266721969F,
      0.9987954562051724F,  0.9996988186962042F,  0.9999247018391445F,
      0.9999811752826011F,  0.9999952938095761F,  0.9999988234517018F,
      0.9999997058628822F,  0.9999999264657178F
  };
  static const float Icoef[14] =
  {   0.0000000000000000F, -1.0000000000000000F, -0.7071067811865474F,
     -0.3826834323650897F, -0.1950903220161282F, -0.0980171403295606F,
     -0.0490676743274180F, -0.0245412285229122F, -0.0122715382857199F,
     -0.0061358846491544F, -0.0030679567629659F, -0.0015339801862847F,
     -0.0007669903187427F, -0.0003834951875714F
  };
  
  nn = N >> 1;
  ie = N;
  for(n=1; n<=LogN; n++) {
    rw = Rcoef[LogN - n];
    iw = Icoef[LogN - n];
    in = ie >> 1;
    ru = 1.0F;
    iu = 0.0F;
    for(j=0; j<in; j++) {
      for(i=j; i<N; i+=ie) {
        io       = i + in;
        rtp      = Rdat[i]  + Rdat[io];
        itp      = Idat[i]  + Idat[io];
        rtq      = Rdat[i]  - Rdat[io];
        itq      = Idat[i]  - Idat[io];
        Rdat[io] = rtq * ru - itq * iu;
        Idat[io] = itq * ru + rtq * iu;
        Rdat[i]  = rtp;
        Idat[i]  = itp;
      }
      sr = ru;
      ru = ru * rw - iu * iw;
      iu = iu * rw + sr * iw;
    }
    ie >>= 1;
  }
  for(j=i=1; i<N; i++) {
    if(i < j) {
      io       = i - 1;
      in       = j - 1;
      rtp      = Rdat[in];
      itp      = Idat[in];
      Rdat[in] = Rdat[io];
      Idat[in] = Idat[io];
      Rdat[io] = rtp;
      Idat[io] = itp;
    }
    k = nn;
    while(k < j){
      j   = j - k;
      k >>= 1;
    }

    j = j + k;
  }
  return true;
}

uint16_t maxFreqAnalysis(float * arr, uint16_t size) {
  float maxValNum = arr[0];
  uint16_t maxValIndex = 0;
  
  for (uint16_t i = 1; i < size; i++) {
    if (arr[i] > maxValNum) {
      maxValNum = arr[i];
      maxValIndex = i;
    }
  }
  return maxValIndex;
}

// UART code
/**@brief   Function for handling app_uart events.
 *
 * @details This function will receive a single character from the app_uart module and append it to 
 *          a string. The string will be be sent over BLE when the last character received was a 
 *          'new line' i.e '\n' (hex 0x0D) or if the string has reached a length of 
 *          @ref NUS_MAX_DATA_LENGTH.
 */
/**@snippet [Handling the data received over UART] */
void uart_error_handle(app_uart_evt_t * p_event){
    if (p_event->evt_type == APP_UART_COMMUNICATION_ERROR)
    {
        APP_ERROR_HANDLER(p_event->data.error_communication);
    }
    else if (p_event->evt_type == APP_UART_FIFO_ERROR)
    {
        APP_ERROR_HANDLER(p_event->data.error_code);
    }
}
/**@snippet [Handling the data received over UART] */


/**@brief  Function for initializing the UART module.
 */
/**@snippet [UART Initialization] */
static void uart_init(void){
    uint32_t                     err_code;
    const app_uart_comm_params_t comm_params =
      {
          RX_PIN_NUMBER,
          TX_PIN_NUMBER,
          RTS_PIN_NUMBER,
          CTS_PIN_NUMBER,
          APP_UART_FLOW_CONTROL_ENABLED,
          false,
          UART_BAUDRATE_BAUDRATE_Baud38400
      };

    APP_UART_FIFO_INIT(&comm_params,
                         UART_RX_BUF_SIZE,
                         UART_TX_BUF_SIZE,
                         uart_error_handle,
                         APP_IRQ_PRIORITY_LOW,
                         err_code);

    APP_ERROR_CHECK(err_code);
}
/**@snippet [UART Initialization] */
      
int main(void) {
	uart_init();
  adcConfig();
  nrf_adc_start();             // trigger next ADC conversion
	//nrf_gpio_cfg_output(24);
	//nrf_gpio_pin_clear(24);
  while (true) {
    __SEV();
    __WFE();
    __WFE();
    if (data1Ready == true || data2Ready == true) {
      //// part to check how FFT works
      /*float  p = 2 * 3.141592653589 / ARRAY_SIZE;              
      for(uint16_t i=0; i<ARRAY_SIZE; i++)
      {
        inArrayR[][i] = sin(p * i);                               // create test signal
        inArrayI[i] = 0.0;                                      
      }*/
      uint8_t arrNum;
      if (data1Ready == true) arrNum = 0;
      else arrNum = 1;
      
      //// FFT and signal analysis
      for(uint16_t i=0; i<ARRAY_SIZE; i++) inArrayI[i] = 0.0;         
      FFT(&inArrayR[arrNum][0], &inArrayI[0], ARRAY_SIZE, LOGN_ARRAY_SIZE);             // calc real and im parts
      static uint16_t fftRes[ARRAY_SIZE/2] = {0};
      for(uint16_t i = 0; i < ARRAY_SIZE/2; i++) {
        fftRes[i] = inArrayR[arrNum][i] = sqrt(inArrayR[arrNum][i]*inArrayR[arrNum][i]+inArrayI[i]*inArrayI[i]);    // calc energy
      }
      
      //// DSP analysis
      static uint16_t maxFreqIndex;
      maxFreqIndex = maxFreqAnalysis(&inArrayR[arrNum][0], ARRAY_SIZE/2);
      
      static uint16_t maxFreq;
      maxFreq = maxFreqIndex * BIN_STEP;
      
      static float maxFreqPower;
      maxFreqPower = inArrayR[arrNum][maxFreqIndex];
        
      static float averPower;
      averPower = 0;
      for(uint16_t i = 0; i < ARRAY_SIZE/2; i++) averPower += inArrayR[arrNum][i];
      averPower = ((averPower - maxFreqPower) / (ARRAY_SIZE/2 - 1));
      
      static uint16_t hapr;
      if (averPower < 0.1) hapr = 0; 
      else hapr = (uint16_t)(((float)maxFreqPower/(float)averPower)*10.0);
      
      // statistics collection
#define STAT_ARRAY_SIZE         256            // size of data array we want to analyze
#define ARRAY_SIZE_2POW         8      
      static uint16_t dataCounter = 0; 
      static uint16_t maxFreqArr[STAT_ARRAY_SIZE], maxFreqPowerArr[STAT_ARRAY_SIZE], 
                      averPowerArr[STAT_ARRAY_SIZE], haprArr[STAT_ARRAY_SIZE], silArr[STAT_ARRAY_SIZE];
      
      maxFreqArr[dataCounter] = maxFreq;
      maxFreqPowerArr[dataCounter] = (uint16_t) maxFreqPower;
      averPowerArr[dataCounter] = (uint16_t) averPower;
      haprArr[dataCounter] = hapr;
      if (averPower < SILENCE_TRSHLD) silArr[dataCounter] = 1;
      else silArr[dataCounter] = 0;
      
      // perform some calculations
      if(++dataCounter > (STAT_ARRAY_SIZE-1)) dataCounter = 0;
      static int32_t maxFreqTt;
      static int32_t maxFreqPowerTt;
      static int32_t averPowerTt;
      static int32_t haprTt;
      static int16_t silTt;
      maxFreqTt = maxFreqPowerTt = averPowerTt = haprTt = silTt = 0;
      
      for (uint16_t i = 0; i < STAT_ARRAY_SIZE; i++) {
        maxFreqTt += maxFreqArr[i];
        maxFreqPowerTt += maxFreqPowerArr[i];
        averPowerTt += averPowerArr[i];
        haprTt += haprArr[i];
        silTt += silArr[i];
      }
      maxFreqTt >>= ARRAY_SIZE_2POW;
      maxFreqPowerTt >>= ARRAY_SIZE_2POW;
      averPowerTt >>= ARRAY_SIZE_2POW;
      haprTt >>= ARRAY_SIZE_2POW;
      
      // let's compare calc results with the pattern
      static uint8_t cryDetectedTimes = 0;
      static bool cryDetected = false;                                 
      if (haprTt > HAPR && averPowerTt > AVER_POWER && (maxFreqTt > LOW_BAND && maxFreqTt < HIGH_BAND) && silTt > NUM_SILENT_SAMPLES) {
        if (cryDetectedTimes < 200) cryDetectedTimes++;
      }
      else {
        if (cryDetectedTimes > 0) cryDetectedTimes--;
      }
      if (cryDetectedTimes > 100){
				cryDetected = true;
				nrf_gpio_pin_clear(24);
			}
      else{
				cryDetected = false;
				nrf_gpio_pin_set(24);
			}
      
      if (arrNum == 0) data1Ready = false;
      else data2Ready = false;
    }
		printf("%f\n", offset);
  }
}

#endif /* NRF51 */

/** @} */
