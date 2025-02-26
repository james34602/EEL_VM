#include "PluginProcessor.h"
#include "PluginEditor.h"
extern "C"
{
#include "simple.h"
#include "../../stb_sprintf.h"
}
void NSEEL_HOSTSTUB_EnterMutex() { }
void NSEEL_HOSTSTUB_LeaveMutex() { }
#include "../../eelCommon.h"
//#include <vld.h>
LiveProgrammableDSP::LiveProgrammableDSP() : treeState(*this, nullptr, "Parameters", createParameterLayout())
{
	fs = 48000.0;
	treeState.addParameterListener("volume", this);
	compileSucessfully = 0;
	codehandleInit = 0;
	codehandleProcess = 0;
	actualLuaVM = 0;
	luaStructRef = LUA_NOREF;
	luaProcessRef = LUA_NOREF;
	NSEEL_start();
	vm = NSEEL_VM_alloc(); // create virtual machine
	vmFs = NSEEL_VM_regvar(vm, "srate");
	*vmFs = fs;
	nSmps = NSEEL_VM_regvar(vm, "nSmps");
	nCh = NSEEL_VM_regvar(vm, "nCh");
	*nCh = JucePlugin_MaxNumInputChannels;
	char splStr[5];
	for (unsigned int i = 0; i < JucePlugin_MaxNumInputChannels; i++)
	{
		sprintf(splStr, "spl%d", i);
		input[i] = NSEEL_VM_regvar(vm, splStr);
	}
	compileContext *ctx = (compileContext*)vm;
	ptrVM_var_name = ctx->varTable_Names;
	ptrVM_var_value = ctx->varTable_Values;
	ptrNumBlocks = &ctx->varTable_numBlocks;
}
void LiveProgrammableDSP::LoadEELCode(char *codeTextInit, char *codeTextProcess, char mode)
{
	ScopedLock lock(criticalSection);
	compileSucessfully = 0;
	compileContext *ctx = (compileContext*)vm;
	NSEEL_VM_freevars(vm);
	NSEEL_init_memRegion(vm);
	memset(ctx->ram_state, 0, sizeof(ctx->ram_state));
	vmFs = NSEEL_VM_regvar(vm, "srate");
	*vmFs = fs;
	nSmps = NSEEL_VM_regvar(vm, "nSmps");
	nCh = NSEEL_VM_regvar(vm, "nCh");
	*nCh = JucePlugin_MaxNumInputChannels;
	char splStr[5];
	for (unsigned int i = 0; i < JucePlugin_MaxNumInputChannels; i++)
	{
		sprintf(splStr, "spl%d", i);
		input[i] = NSEEL_VM_regvar(vm, splStr);
	}
	ptrVM_var_name = ctx->varTable_Names;
	ptrVM_var_value = ctx->varTable_Values;
	ptrNumBlocks = &ctx->varTable_numBlocks;
	if (codehandleInit)
	{
		NSEEL_code_free(codehandleInit);
		codehandleInit = 0;
	}
	if (codehandleProcess)
	{
		NSEEL_code_free(codehandleProcess);
		codehandleProcess = 0;
	}
	ctx->functions_common = 0;
	codehandleInit = NSEEL_code_compile_ex(vm, codeTextInit, 0, 1);
	EEL_STRING_STDOUT_WRITE(ctx->last_error_string, strlen(ctx->last_error_string));
	if (!codehandleInit)
	{
		compileSucessfully = 0;
		return;
	}
	NSEEL_code_execute(codehandleInit);
	codehandleProcess = NSEEL_code_compile(vm, codeTextProcess, 0);
	EEL_STRING_STDOUT_WRITE(ctx->last_error_string, strlen(ctx->last_error_string));
	if (codehandleInit && codehandleProcess)
	{
		if (codehandleProcess)
			compileSucessfully = mode;
		else
			compileSucessfully = 0;
	}
	else
		compileSucessfully = 0;
	char *msg = "EEL: Code looks OK\nLua compilation will not be run\n";
	EEL_STRING_STDOUT_WRITE(msg, strlen(msg));
}
LiveProgrammableDSP::~LiveProgrammableDSP()
{
	if (vm)
	{
		NSEEL_code_free(codehandleInit);
		NSEEL_code_free(codehandleProcess);
		NSEEL_VM_free(vm);
	}
	NSEEL_quit();
	if (actualLuaVM)
	{
		if (luaProcessRef != LUA_NOREF)
			luaL_unref(actualLuaVM, LUA_REGISTRYINDEX, luaProcessRef);
		if (luaStructRef != LUA_NOREF)
			luaL_unref(actualLuaVM, LUA_REGISTRYINDEX, luaStructRef);
		lua_close(actualLuaVM);
	}
}
AudioProcessorValueTreeState::ParameterLayout LiveProgrammableDSP::createParameterLayout()
{
	std::vector<std::unique_ptr<RangedAudioParameter>> parameters;
	parameters.push_back(std::make_unique<AudioParameterFloat>("volume", "volume", NormalisableRange<float>(0.0, 1.0, 0.01), 1.0, "val"));
	return { parameters.begin(), parameters.end() };
}
void LiveProgrammableDSP::parameterChanged(const String& parameter, float newValue)
{
	volume = newValue;
}
void LiveProgrammableDSP::prepareToPlay(double sampleRate, int samplesPerBlock)
{
	fs = (double)sampleRate;
	*vmFs = fs;
	if (compileSucessfully == 1 || compileSucessfully == 2)
	{
		ScopedLock lock(criticalSection);
		NSEEL_code_execute(codehandleInit);
	}
	if (compileSucessfully == 3)
	{
		ScopedLock lock(criticalSection);
		int top = lua_gettop(actualLuaVM);
		int retType = lua_getglobal(actualLuaVM, "init");  // function to be called
		lua_pushnumber(actualLuaVM, fs);   // push 1st argument
		// do the call (1 arguments, 1 table)
		lua_call(actualLuaVM, 1, LUA_MULTRET);
		// check the return value
		if (lua_gettop(actualLuaVM) - top) {
			// store userdata to a pointer
			luaStructRef = luaL_ref(actualLuaVM, LUA_REGISTRYINDEX);
		}
		else
		{
			char erMsg[70] = "Lua: Fail to get reference while re-set sample rate, stop compiling\n";
			EEL_STRING_STDOUT_WRITE(erMsg, strlen(erMsg));
			lua_close(actualLuaVM);
		}
	}
}
const String LiveProgrammableDSP::getName() const
{
	return JucePlugin_Name;
}
bool LiveProgrammableDSP::acceptsMidi() const
{
#if JucePlugin_WantsMidiInput
	return true;
#else
	return false;
#endif
}
bool LiveProgrammableDSP::producesMidi() const
{
#if JucePlugin_ProducesMidiOutput
	return true;
#else
	return false;
#endif
}
bool LiveProgrammableDSP::isMidiEffect() const
{
#if JucePlugin_IsMidiEffect
	return true;
#else
	return false;
#endif
}
double LiveProgrammableDSP::getTailLengthSeconds() const
{
	return 0.0;
}
int LiveProgrammableDSP::getNumPrograms()
{
	return 1;
}
int LiveProgrammableDSP::getCurrentProgram()
{
	return 0;
}
void LiveProgrammableDSP::setCurrentProgram(int)
{
}
const String LiveProgrammableDSP::getProgramName(int)
{
	return {};
}
void LiveProgrammableDSP::changeProgramName(int, const String&)
{
}
void LiveProgrammableDSP::releaseResources()
{
}
// ========================= Define processBlock Method =========================
void LiveProgrammableDSP::processBlock(AudioBuffer<float>& buffer, MidiBuffer&)
{
	// number of samples per buffer
	const int n = buffer.getNumSamples();
	const float *inputs[JucePlugin_MaxNumInputChannels];
	float* outputs[JucePlugin_MaxNumInputChannels];
	for (int i = 0; i < JucePlugin_MaxNumInputChannels; i++)
	{
		inputs[i] = buffer.getReadPointer(i);
		outputs[i] = buffer.getWritePointer(i);
	}
	if (compileSucessfully == 1)
	{
		ScopedLock lock(criticalSection);
		*nSmps = 1;
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < JucePlugin_MaxNumInputChannels; j++)
				*input[j] = inputs[j][i];
			NSEEL_code_execute(codehandleProcess);
			for (int j = 0; j < JucePlugin_MaxNumInputChannels; j++)
				outputs[j][i] = *input[j];
		}
	}
	if (compileSucessfully == 2)
	{
		*nSmps = n;
		float *dat = dataSectionToRamDisk(vm, JucePlugin_MaxNumInputChannels * n);
		for (int j = 0; j < JucePlugin_MaxNumInputChannels; j++)
			for (int i = 0; i < n; i++)
				dat[i + n * j] = inputs[j][i];
		{
			ScopedLock lock(criticalSection);
			NSEEL_code_execute(codehandleProcess);
		}
		for (int j = 0; j < JucePlugin_MaxNumInputChannels; j++)
			for (int i = 0; i < n; i++)
				outputs[j][i] = dat[i + n * j];
	}
	if (compileSucessfully == 3)
	{
		ScopedLock lock(criticalSection);
		for (int i = 0; i < n; i++)
		{
			lua_rawgeti(actualLuaVM, LUA_REGISTRYINDEX, luaProcessRef);
			lua_rawgeti(actualLuaVM, LUA_REGISTRYINDEX, luaStructRef);
			// Push variable to stack
			lua_pushnumber(actualLuaVM, inputs[0][i]);
			lua_pushnumber(actualLuaVM, inputs[1][i]);
			lua_pushnumber(actualLuaVM, inputs[2][i]);
			lua_pushnumber(actualLuaVM, inputs[3][i]);
			lua_pushnumber(actualLuaVM, inputs[4][i]);
			lua_pushnumber(actualLuaVM, inputs[5][i]);
			lua_call(actualLuaVM, 7, 6);
			outputs[0][i] = (float)lua_tonumber(actualLuaVM, -6);
			outputs[1][i] = (float)lua_tonumber(actualLuaVM, -5);
			outputs[2][i] = (float)lua_tonumber(actualLuaVM, -4);
			outputs[3][i] = (float)lua_tonumber(actualLuaVM, -3);
			outputs[4][i] = (float)lua_tonumber(actualLuaVM, -2);
			outputs[5][i] = (float)lua_tonumber(actualLuaVM, -1);
			lua_pop(actualLuaVM, 6);
		}
	}
}
bool LiveProgrammableDSP::hasEditor() const
{
	return true;
}
AudioProcessorEditor* LiveProgrammableDSP::createEditor()
{
	return new LiveProgEditor(*this);
}
void LiveProgrammableDSP::getStateInformation(MemoryBlock& destData)
{
	MemoryOutputStream stream(destData, true);
	unsigned char *compressed;
	size_t sz = 0;
	const char *sampleData = restoredText.toRawUTF8();
	size_t inLen = strlen(sampleData);
	if (!inLen)
		return;
	int rc = simpleCompress(ELZMA_lzma, (unsigned char *)sampleData, inLen, &compressed, &sz);
	if (rc != ELZMA_E_OK)
		return;
	stream.write(compressed, sz);
	free(compressed);
//	stream.writeString(restoredText);
}
void LiveProgrammableDSP::setStateInformation(const void* data, int sizeInBytes)
{
	MemoryInputStream stream(data, static_cast<size_t> (sizeInBytes), false);
	unsigned char *compressed = (unsigned char*)malloc(sizeInBytes);
	memset(compressed, 0, sizeInBytes);
	stream.read((void*)compressed, sizeInBytes);
	size_t sz = 0;
	unsigned char *decompressed;
	int rc = simpleDecompress(ELZMA_lzma, compressed, sizeInBytes, &decompressed, &sz);
	free(compressed);
	if (rc != ELZMA_E_OK)
		return;
	if (!sz)
	{
		free(decompressed);
		return;
	}
	decompressed = (unsigned char*)realloc(decompressed, sz + 1);
	decompressed[sz] = 0;
	restoredText = String(CharPointer_UTF8((char*)decompressed));
	free(decompressed);
	playButtonClicked(restoredText);
}
int LiveProgrammableDSP::compileEEL(const char *eelCode, size_t strLen)
{
	const char *initSegment = strstr(eelCode, "@init");
	if (!initSegment)
	{
		char *errorMsg = "EEL: @init section not found\n";
		EEL_STRING_STDOUT_WRITE(errorMsg, 29);
		return 0;
	}
	else
		initSegment += 6;
	char mode = 0;
	const char *processSegment = strstr(eelCode, "@frame");
	if (processSegment)
	{
		processSegment += 7;
		mode = 0;
	}
	else
	{
		processSegment = strstr(eelCode, "@sample");
		if (!processSegment)
		{
			char *errorMsg = "EEL: @sample section not found\n";
			EEL_STRING_STDOUT_WRITE(errorMsg, 31);
			return 0;
		}
		else
		{
			processSegment += 8;
			mode = 1;
		}
	}
	char *codeTextInit = (char*)malloc(strLen * sizeof(char));
	char *codeTextProcess = (char*)malloc(strLen * sizeof(char));
	memset(codeTextInit, 0, strLen * sizeof(char));
	memset(codeTextProcess, 0, strLen * sizeof(char));
	if (initSegment && processSegment)
	{
		if (initSegment < processSegment)
		{
			const char *sampleSegment = strstr(eelCode, "@sample");
			if ((sampleSegment - processSegment) > 0)
			{
				int cpyLen = processSegment - initSegment - 8;
				if (cpyLen > 0)
					strncpy(codeTextInit, initSegment, cpyLen);
				if (processSegment - eelCode < strLen)
					strncpy(codeTextProcess, processSegment, (sampleSegment - processSegment) - 1);
			}
			else
			{
				int cpyLen;
				if (sampleSegment)
					cpyLen = sampleSegment - initSegment;
				else
					cpyLen = processSegment - initSegment - (7 + 1);
				if (cpyLen > 0)
					strncpy(codeTextInit, initSegment, cpyLen);
				if (processSegment - eelCode < strLen)
					strcpy(codeTextProcess, processSegment);
			}
			LoadEELCode(codeTextInit, codeTextProcess, mode);
		}
		else
		{
			free(codeTextInit);
			free(codeTextProcess);
			return 0;
		}
	}
	free(codeTextInit);
	free(codeTextProcess);
	return 1;
}
int LiveProgrammableDSP::validateLua(const char *luaCode, size_t strLen)
{
	lua_State *L = luaL_newstate();  // create state
	luaL_openlibs(L);  // open standard libraries
	lua_settop(L, 0);
	char *concatedStr = (char*)malloc(strLen + 256);
	int outLen = stbsp_snprintf(concatedStr, strLen + 256, "%s\nprocess(init(%lf), %lf, %lf)", luaCode, nseel_int_rand(192000.0), nseel_int_rand(1.0), nseel_int_rand(1.0));
	int err = luaL_loadbuffer(L, concatedStr, outLen, "LiveProg");
	free(concatedStr);
	EEL_STRING_STDOUT_WRITE("Lua: Test run\n", 14);
	if (err)
	{
		char erMsg[256] = { 0 };
		int msgLen = stbsp_snprintf(erMsg, 256, "Lua error %d\n%s\n", err, lua_tostring(L, -1));
		EEL_STRING_STDOUT_WRITE(erMsg, msgLen);
		lua_pop(L, 1);  // pop error message from the stack
		lua_close(L);
		return -1;
	}
	err = lua_pcall(L, 0, 0, 0);
	if (err)
	{
		char erMsg[256] = { 0 };
		int msgLen = stbsp_snprintf(erMsg, 256, "Lua error %d\n%s\n", err, lua_tostring(L, -1));
		EEL_STRING_STDOUT_WRITE(erMsg, msgLen);
		lua_pop(L, 1);  // pop error message from the stack
		lua_close(L);
		return -1;
	}
	// Table return
	int userDataRef = LUA_NOREF;
	int top = lua_gettop(L);
	int retType = lua_getglobal(L, "init");  // function to be called
	lua_pushnumber(L, fs);   // push 1st argument
	// do the call (1 arguments, 1 table)
	if (lua_pcall(L, 1, LUA_MULTRET, 0) != 0)
	{
		char erMsg[256] = { 0 };
		int msgLen = stbsp_snprintf(erMsg, 256, "Lua error %d\nerror running function 'init': %s\n", err, lua_tostring(L, -1));
		EEL_STRING_STDOUT_WRITE(erMsg, msgLen);
		lua_pop(L, 1);  // pop error message from the stack
		lua_close(L);
		return 1;
	}
	int retClass = lua_istable(L, -1);
	if (!retClass)
	{
		char erMsg[55] = "Lua: init function must return table, stop compiling\n";
		EEL_STRING_STDOUT_WRITE(erMsg, 54);
		lua_close(L);
		return 2;
	}
	// check the return value
	if (lua_gettop(L) - top) {
		// store userdata to a pointer
		userDataRef = luaL_ref(L, LUA_REGISTRYINDEX);
	}
	else
	{
		char erMsg[45] = "Lua: Fail to get reference, stop compiling\n";
		EEL_STRING_STDOUT_WRITE(erMsg, 44);
		lua_close(L);
		return 3;
	}
	lua_getglobal(L, "process");  // function to be called
	int process = luaL_ref(L, LUA_REGISTRYINDEX);
	if (userDataRef != LUA_NOREF && userDataRef != LUA_REFNIL)
	{
		lua_rawgeti(L, LUA_REGISTRYINDEX, process);
		lua_rawgeti(L, LUA_REGISTRYINDEX, userDataRef);
		// Push variable to stack
		lua_pushnumber(L, 0.5);
		lua_pushnumber(L, -0.15);
		if (lua_pcall(L, 3, 2, 0) != 0)
		{
			char erMsg[256] = { 0 };
			int msgLen = stbsp_snprintf(erMsg, 256, "Lua error %d\nerror running function 'process': %s\n", err, lua_tostring(L, -1));
			EEL_STRING_STDOUT_WRITE(erMsg, msgLen);
			lua_pop(L, 1);  // pop error message from the stack
			lua_close(L);
			return 1;
		}
		double y1 = lua_tonumber(L, -2);
		double y2 = lua_tonumber(L, -1);
		lua_pop(L, 2);
	}
	else
	{
		char erMsg[42] = "Lua: Reference go wrong, stop compiling\n";
		EEL_STRING_STDOUT_WRITE(erMsg, 41);
	}
	luaL_unref(L, LUA_REGISTRYINDEX, process);
	luaL_unref(L, LUA_REGISTRYINDEX, userDataRef);
	lua_close(L);
	return LUA_OK;
}
void LiveProgrammableDSP::LoadLuaCode(const char *luaCode, size_t strLen)
{
	ScopedLock lock(criticalSection);
	compileSucessfully = 0;
	int i;
	if (actualLuaVM)
	{
		if (luaProcessRef != LUA_NOREF)
			luaL_unref(actualLuaVM, LUA_REGISTRYINDEX, luaProcessRef);
		luaProcessRef = LUA_NOREF;
		if (luaStructRef != LUA_NOREF)
			luaL_unref(actualLuaVM, LUA_REGISTRYINDEX, luaStructRef);
		luaStructRef = LUA_NOREF;
		lua_close(actualLuaVM);
		actualLuaVM = 0;
	}
	actualLuaVM = luaL_newstate();  // create state
	luaL_openlibs(actualLuaVM);  // open standard libraries
	lua_settop(actualLuaVM, 0);
	int err = luaL_loadbuffer(actualLuaVM, luaCode, strLen, "LiveProg");
	if (err)
	{
		char erMsg[256] = { 0 };
		int msgLen = stbsp_snprintf(erMsg, 256, "Lua error %d\n%s\nFatal code loading error after validation\n", err, lua_tostring(actualLuaVM, -1));
		EEL_STRING_STDOUT_WRITE(erMsg, msgLen);
		lua_pop(actualLuaVM, 1);  // pop error message from the stack
		lua_close(actualLuaVM);
		return;
	}
	lua_call(actualLuaVM, 0, 0);
	int top = lua_gettop(actualLuaVM);
	int retType = lua_getglobal(actualLuaVM, "init");  // function to be called
	lua_pushnumber(actualLuaVM, fs);   // push 1st argument
	// do the call (1 arguments, 1 table)
	lua_call(actualLuaVM, 1, LUA_MULTRET);
	// check the return value
	if (lua_gettop(actualLuaVM) - top) {
		// store userdata to a pointer
		luaStructRef = luaL_ref(actualLuaVM, LUA_REGISTRYINDEX);
	}
	else
	{
		char erMsg[105] = "Lua: Fail to get reference after validation, stop compiling\nFatal get reference error after validation\n";
		EEL_STRING_STDOUT_WRITE(erMsg, 104);
		lua_close(actualLuaVM);
		return;
	}
	lua_getglobal(actualLuaVM, "process");  // function to be called
	luaProcessRef = luaL_ref(actualLuaVM, LUA_REGISTRYINDEX);
	compileSucessfully = 3;
	EEL_STRING_STDOUT_WRITE("Lua: Code looks OK\n", 16);
}
int LiveProgrammableDSP::compileLua(const char *luaCode, size_t strLen)
{
	EEL_STRING_STDOUT_WRITE("Code text is not EEL, compiling code as Lua\n", 45);
	int ret = validateLua(luaCode, strLen);
	if (ret != LUA_OK)
	{
		EEL_STRING_STDOUT_WRITE("Lua compilation failed\n", 24);
		compileSucessfully = 0;
		return 0;
	}
	else
	{
		LoadLuaCode(luaCode, strLen);
		return 1;
	}
}
void LiveProgrammableDSP::playButtonClicked(String codetext)
{
	int strLen = codetext.length();
	const char *textbuf = codetext.toRawUTF8();
	int eelSuccess = compileEEL(textbuf, strLen);
	if (!eelSuccess)
		int luaSuccess = compileLua(textbuf, strLen);
}
AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
	return new LiveProgrammableDSP();
}
