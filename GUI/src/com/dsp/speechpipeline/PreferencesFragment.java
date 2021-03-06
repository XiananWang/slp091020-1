package com.dsp.speechpipeline;

import java.util.Locale;

import android.content.SharedPreferences;
import android.content.SharedPreferences.OnSharedPreferenceChangeListener;
import android.os.Bundle;
import android.preference.ListPreference;
import android.preference.Preference;
import android.preference.Preference.OnPreferenceChangeListener;
import android.preference.PreferenceFragment;

public class PreferencesFragment extends PreferenceFragment implements OnSharedPreferenceChangeListener {
	
	private Monitor main;

	public void onCreate(Bundle savedInstanceState) {
		super.onCreate(savedInstanceState);
		
		main = Settings.getCallbackInterface();
		
		//preference options are defined in res/xml/prefs.xml
		//preference names are defined in res/values/strings.xml
 		addPreferencesFromResource(R.xml.prefs);
 		
 		//sound output option
 		Preference preference = findPreference(getString(R.string.prefPlayback));
 		preference.setSummary("Current: " + (Settings.playback?"Enabled":"Disabled"));
 		
 		//sound channel option
 		preference = findPreference(getString(R.string.prefOutputStream));
 		preference.setSummary("Current: " + Settings.getOutput());
 		
 		//sampling frequency options
        preference = findPreference(getString(R.string.prefSamplingFreq));
        ((ListPreference)preference).setEntries(Settings.samplingRates);
        ((ListPreference)preference).setEntryValues(Settings.samplingRateValues);
        preference.setSummary("Current: " + Settings.Fs + "Hz");
        
        //window time option
        preference = findPreference(getString(R.string.prefWindowTime));
        preference.setSummary("Current: " + Settings.windowTime + "ms");
        preference.setOnPreferenceChangeListener(preferenceChange);
        
        //step time option
        preference = findPreference(getString(R.string.prefStepTime));
        preference.setSummary("Current: " + Settings.stepTime + "ms");
        preference.setOnPreferenceChangeListener(preferenceChange);
        
        //debug output option
        preference = findPreference(getString(R.string.prefDebug));
        preference.setSummary("Current: " + Settings.getDebugLevel());
        
        //debug text output option
        preference = findPreference(getString(R.string.prefDebugOutput));
        preference.setSummary("Current: " + Settings.getDebugOutput());
	}

	public void onResume() {
		super.onResume();
		getPreferenceManager().getSharedPreferences().registerOnSharedPreferenceChangeListener(this);
	}
	
	public void onPause() {
		getPreferenceManager().getSharedPreferences().unregisterOnSharedPreferenceChangeListener(this);
		super.onPause();
	}
	
	private OnPreferenceChangeListener preferenceChange = new OnPreferenceChangeListener() {
		public boolean onPreferenceChange(Preference preference, Object changedValue) {
			String testString = changedValue.toString().trim();
			if(testString == "") {
				return false;
			} else {
				try {
					Float.parseFloat(testString);
					return true;
				} catch(NumberFormatException e){}
				return false;
			}
		}
	};

	public void onSharedPreferenceChanged(SharedPreferences sharedPreferences, String key) {
		
		Preference preference = findPreference(key);
		
		//sound output option
		if (key.equals(getString(R.string.prefPlayback))) {
			boolean setting = sharedPreferences.getBoolean(key, false);
			if(Settings.setPlayback(setting)) {
				String[] result = {"Playback disabled.", "Playback enabled for " + Settings.getOutput().toLowerCase(Locale.US) + " output."};
				preference.setSummary("Current: " + (Settings.playback?"Enabled":"Disabled"));
				preference = findPreference(getString(R.string.prefOutputStream));
		 		preference.setSummary("Current: " + Settings.getOutput()); 	
		 		main.notify(result[Settings.playback?1:0]);
			}
		}
		
		//sound channel option
		else if (key.equals(getString(R.string.prefOutputStream))) {
			int setting = Integer.parseInt(sharedPreferences.getString(key, "2"));
			if(Settings.setOutput(setting)) {
				preference.setSummary("Current: " + Settings.getOutput());
				String[] result = {"Playback set to original output.", "Playback set to filtered output."};
				main.notify(result[Settings.output?1:0]);
			}
		}
		
		//sampling frequency options
		else if (key.equals(getString(R.string.prefSamplingFreq))) {
			int setting = Integer.parseInt(sharedPreferences.getString(key, "8000"));
			if(Settings.setSamplingFrequency(setting)) {
				preference.setSummary("Current: " + Settings.Fs + "Hz");
				main.notify("Sampling rate set to " + Settings.Fs + "Hz.");
			}
		}
		
		//window time option
		else if (key.equals(getString(R.string.prefWindowTime))) {
			float setting = Float.parseFloat(sharedPreferences.getString(key, "11.0"));
			if(Settings.setWindowSize(setting)) {
				preference.setSummary("Current: " + Settings.windowTime + "ms");
				main.notify("Window time set to " + Settings.windowTime + "ms.");
			}
		}
		
		//step time option
		else if (key.equals(getString(R.string.prefStepTime))) {
			float setting = Float.parseFloat(sharedPreferences.getString(key, "5.0"));
			if(Settings.setStepSize(setting)) {
				preference.setSummary("Current: " + Settings.stepTime + "ms");
				main.notify("Step time set to " + Settings.stepTime + "ms.");
			}
		}
		
		//debug output option
		else if (key.equals(getString(R.string.prefDebug))) {
			int setting = Integer.parseInt(sharedPreferences.getString(key, "4"));
			if(Settings.setDebugLevel(setting)) {
				preference.setSummary("Current: " + Settings.getDebugLevel());
				preference = findPreference(getString(R.string.prefDebugOutput));
		 		preference.setSummary("Current: " + Settings.getDebugOutput()); 	
				String[] result = {"All debug outputs enabled.", "Text file output enabled.", "PCM ouput enabled.", "Debug ouput disabled."};
				main.notify(result[Settings.debugLevel-1]);
			}
		}
		
		//debug text output option
		else if (key.equals(getString(R.string.prefDebugOutput))) {
			int setting = Integer.parseInt(sharedPreferences.getString(key, "0"));
			if(Settings.setDebugOutput(setting)) {
				preference.setSummary("Current: " + Settings.getDebugOutput());
				main.notify( Settings.getDebugOutput() + " text output selected.");
			}
		}
	}
}