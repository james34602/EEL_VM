#pragma once
#include "../JuceLibraryCode/JuceHeader.h"
#include "../../eelCommon.h"
extern "C"
{
#include "lua/lua.h"
#include "lua/lauxlib.h"
#include "lua/lualib.h"
}
class LiveProgrammableDSP : public AudioProcessor, public AudioProcessorValueTreeState::Listener
{
public:
	// declare constractor and deconstractor
	LiveProgrammableDSP();
	~LiveProgrammableDSP();
	// declare default methods
	void prepareToPlay(double sampleRate, int samplesPerBlock) override;
	void releaseResources() override;
	void parameterChanged(const String& parameter, float newValue) override;
	void processBlock(AudioBuffer<float> &buffer, MidiBuffer &) override;

	// declare extra default methods
	AudioProcessorEditor* createEditor() override;
	bool hasEditor() const override;
	const String getName() const override;
	bool acceptsMidi() const override;
	bool producesMidi() const override;
	bool isMidiEffect() const override;
	double getTailLengthSeconds() const override;
	int getNumPrograms() override;
	int getCurrentProgram() override;
	void setCurrentProgram(int index) override;
	const String getProgramName(int index) override;
	void changeProgramName(int index, const String& newName) override;
	void getStateInformation(MemoryBlock& destData) override;
	void setStateInformation(const void* data, int sizeInBytes) override;
	// declare slider parameter vars
	AudioProcessorValueTreeState treeState;
	AudioProcessorValueTreeState::ParameterLayout createParameterLayout();
	// DSP variable
	NSEEL_VMCTX vm;
	NSEEL_CODEHANDLE codehandleInit, codehandleProcess;
	float *vmFs, *input[JucePlugin_MaxNumInputChannels], *nSmps, *nCh;
	// Lua VM
	lua_State *actualLuaVM;
	int luaStructRef, luaProcessRef;
	//
	void LoadEELCode(char *codeTextInit, char *codeTextProcess, char mode);
	void LoadLuaCode(const char *eelCode, size_t strLen);
	int compileEEL(const char *eelCode, size_t strLen);
	int validateLua(const char *eelCode, size_t strLen);
	int compileLua(const char *eelCode, size_t strLen);
	void playButtonClicked(String codetext);
	String restoredText;

	char ***ptrVM_var_name;
	float **ptrVM_var_value;
	int *ptrNumBlocks;
private:
	char *codeTextInit, *codeTextProcess;
	CriticalSection criticalSection;
	float volume;
	int compileSucessfully;
	double fs;
	// Juce extra detection
	JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(LiveProgrammableDSP)
};