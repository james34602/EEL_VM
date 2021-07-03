-- Meta class
BiquadSOS =
{ b0 = 0, b1 = 0, b2 = 0, a1 = 0, a2 = 0,
  v1L = 0, v2L = 0, v1R = 0, v2R = 0 }
-- Derived class method new
function BiquadSOS:new(o, coeffB0, coeffB1, coeffB2, coeffA1, coeffA2)
	o = o or {}
	setmetatable(o, self)
	self.__index = self
	self.b0 = coeffB0
	self.b1 = coeffB1
	self.b2 = coeffB2
	self.a1 = coeffA1
	self.a2 = coeffA2
	self:clearState()
	return o
end
-- Derived class method printCoefficients
function BiquadSOS:printCoefficients()
	io.write("fvtool([", self.b0, ",", self.b1, ",", self.b2, "],[1,", self.a1, ",", self.a2, "]);")
end
-- Derived class method clearState
function BiquadSOS:clearState()
	self.v1L = 0
	self.v2L = 0
	self.v1R = 0
	self.v2R = 0
end
-- Lua version of https://github.com/james34602/JamesDSPManager/blob/master/Open_source_edition/Audio_Engine/eclipse_libjamesdsp_free_bp/jni/vdc.c
function BiquadSOS:Process(x1, x2)
	local w1 = x1 - self.a1 * self.v1L - self.a2 * self.v2L;
	local y1 = self.b0 * w1 + self.b1 * self.v1L + self.b2 * self.v2L;
	self.v2L = self.v1L;
	self.v1L = w1;
	local w2 = x2 - self.a1 * self.v1R - self.a2 * self.v2R;
	local y2 = self.b0 * w2 + self.b1 * self.v1R + self.b2 * self.v2R;
	self.v2R = self.v1R;
	self.v1R = w2;
	return y1, y2
end
-- Class end

function designPeakingFilter(dbGain, centreFreq, fs, dBandwidthOrQOrS)
	if (centreFreq <= math.eps or fs <= math.eps) then
		return 0, 0, 0, 0, 0
	end
	local ln2div2 = math.log(2.0) / 2.0;
	local lingain = math.pow(10.0, dbGain / 40.0);
	local omega = (6.2831853071795862 * centreFreq) / fs;
	local num3 = math.sin(omega);
	local cs = math.cos(omega);
	local alpha = num3 * math.sinh((ln2div2 * dBandwidthOrQOrS * omega) / num3);
	local B0 = 1.0 + (alpha * lingain);
	local B1 = -2.0 * cs;
	local B2 = 1.0 - (alpha * lingain);
	local A0 = 1.0 + (alpha / lingain);
	local A1 = -2.0 * cs;
	local A2 = 1.0 - (alpha / lingain);
	return (B0 / A0), (B1 / A0), (B2 / A0), (A1 / A0), (A2 / A0)
end

function init(srate) -- Test passed, ***Reference this function's own return, a table return!!!
	local b0, b1, b2, a1, a2 = designPeakingFilter(5, 1000, srate, 2);
	local myClassIIR = BiquadSOS:new(nil, b0, b1, b2, a1, a2)
--	io.write("fvtool([", b0, ",", b1, ",", b2, "],[1,", a1, ",", a2, "]);")
	myClassIIR:printCoefficients()
	-- Construct some members
	local mTbl = {
	gainTest = 0.41,
	claszObj = myClassIIR
	}
	return mTbl
end

function process(CReferenceTbl, spl0, spl1, spl2, spl3, spl4, spl5) -- Passing the reference from C to another Lua function
	spl0, spl1 = CReferenceTbl.claszObj:Process(spl0, spl1)
	return spl0, spl1, spl2, spl3, spl4, spl5
end