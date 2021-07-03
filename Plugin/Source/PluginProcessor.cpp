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
	input1 = NSEEL_VM_regvar(vm, "spl0");
	input2 = NSEEL_VM_regvar(vm, "spl1");
	input3 = NSEEL_VM_regvar(vm, "spl2");
	input4 = NSEEL_VM_regvar(vm, "spl3");
	input5 = NSEEL_VM_regvar(vm, "spl4");
	input6 = NSEEL_VM_regvar(vm, "spl5");
	compileContext *ctx = (compileContext*)vm;
	ptrVM_var_name = ctx->varTable_Names;
	ptrVM_var_value = ctx->varTable_Values;
	ptrNumBlocks = &ctx->varTable_numBlocks;
}
void LiveProgrammableDSP::LoadEELCode(char *codeTextInit, char *codeTextProcess)
{
	ScopedLock lock(criticalSection);
	compileSucessfully = 0;
	compileContext *ctx = (compileContext*)vm;
	NSEEL_VM_freevars(vm);
	NSEEL_init_string(vm);
	vmFs = NSEEL_VM_regvar(vm, "srate");
	*vmFs = fs;
	input1 = NSEEL_VM_regvar(vm, "spl0");
	input2 = NSEEL_VM_regvar(vm, "spl1");
	input3 = NSEEL_VM_regvar(vm, "spl2");
	input4 = NSEEL_VM_regvar(vm, "spl3");
	input5 = NSEEL_VM_regvar(vm, "spl4");
	input6 = NSEEL_VM_regvar(vm, "spl5");
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
			compileSucessfully = 1;
		else
			compileSucessfully = 0;
	}
	else
	{
		compileSucessfully = 0;
	}
	char *msg = "EEL: Code looks OK\nLua compilation will not be runned\n";
	EEL_STRING_STDOUT_WRITE(msg, 56);
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
	if (compileSucessfully == 1)
	{
		ScopedLock lock(criticalSection);
		NSEEL_code_execute(codehandleInit);
	}
	if (compileSucessfully == 2)
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
			EEL_STRING_STDOUT_WRITE(erMsg, 69);
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
	// input channels
	const float *inputs[6] = { buffer.getReadPointer(0), buffer.getReadPointer(1), buffer.getReadPointer(2), buffer.getReadPointer(3), buffer.getReadPointer(4), buffer.getReadPointer(5) };
	// output channels
	float* const outputs[6] = { buffer.getWritePointer(0), buffer.getWritePointer(1), buffer.getWritePointer(2), buffer.getWritePointer(3), buffer.getWritePointer(4), buffer.getWritePointer(5) };
	if (compileSucessfully == 1)
	{
		ScopedLock lock(criticalSection);
		for (int i = 0; i < n; i++)
		{
			*input1 = inputs[0][i];
			*input2 = inputs[1][i];
			*input3 = inputs[2][i];
			*input4 = inputs[3][i];
			*input5 = inputs[4][i];
			*input6 = inputs[5][i];
			NSEEL_code_execute(codehandleProcess);
			outputs[0][i] = (float)*input1;
			outputs[1][i] = (float)*input2;
			outputs[2][i] = (float)*input3;
			outputs[3][i] = (float)*input4;
			outputs[4][i] = (float)*input5;
			outputs[5][i] = (float)*input6;
		}
	}
	if (compileSucessfully == 2)
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
//	restoredText = stream.readString();
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
	const char *processSegment = strstr(eelCode, "@sample");
	if (!processSegment)
	{
		char *errorMsg = "EEL: @sample section not found\n";
		EEL_STRING_STDOUT_WRITE(errorMsg, 31);
		return 0;
	}
	else
		processSegment += 8;
	char *codeTextInit = (char*)malloc(strLen * sizeof(char));
	char *codeTextProcess = (char*)malloc(strLen * sizeof(char));
	memset(codeTextInit, 0, strLen * sizeof(char));
	memset(codeTextProcess, 0, strLen * sizeof(char));
	if (initSegment && processSegment)
	{
		if (initSegment < processSegment)
		{
			int cpyLen = processSegment - initSegment - (8 + 1);
			if (cpyLen > 0)
				strncpy(codeTextInit, initSegment, cpyLen);
			if (processSegment - eelCode < strLen)
			{
				strcpy(codeTextProcess, processSegment);
			}
			LoadEELCode(codeTextInit, codeTextProcess);
		}
		else
		{
			strcpy(codeTextInit, initSegment);
			int cpyLen = initSegment - processSegment - (6 + 1);
			if (cpyLen > 0)
				strncpy(codeTextProcess, processSegment, cpyLen);
			LoadEELCode(codeTextInit, codeTextProcess);
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
	compileSucessfully = 2;
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
