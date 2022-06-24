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
void LiveProgrammableDSP::playButtonClicked(String codetext)
{
	int strLen = codetext.length();
	const char *textbuf = codetext.toRawUTF8();
	int eelSuccess = compileEEL(textbuf, strLen);
}
AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
	return new LiveProgrammableDSP();
}
