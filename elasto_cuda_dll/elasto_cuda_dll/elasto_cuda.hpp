
#ifndef   ELASTO_CUDA

#define   ELASTO_CUDA

#ifndef      ELASTO__EXPORTS

#define     ELASTO_API __declspec(dllexport) //导出符号宏定义


#else

#define     ELASTO_API __declspec(dllimport) //导出符号宏定义

#endif

#include  <string>


 // #include  "lib_add_cuda.cuh"

// #include  "CElasto.h"  

//#include  "SysConfig.h"


struct  EInput    ;

struct  EOutput   ;

struct ConfigParam;

class   CudaMain ;


const  char DefaultElastoConfFile[] = ".\\config.ini";


ELASTO_API    class   ElastoCuda    {

	
private :

	std::string  *  initFile;

	CudaMain  *    cudaMain;

	ConfigParam *  config;



private:                         //读取原始RF数据 ，读取滤波器数据，读取配置数据


	void   readRFData(const EInput & in);








public:

	ElastoCuda();



	~ElastoCuda(); 

	virtual bool  isAvailable();

	virtual int   ci() const { return 40 };

	virtual void  init(const std::string &ini_file);

	virtual bool  process( const EInput &input, EOutput &output);










};












#endif