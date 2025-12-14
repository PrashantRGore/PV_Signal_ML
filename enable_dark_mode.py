from pathlib import Path

# Create .streamlit config directory
config_dir = Path('.streamlit')
config_dir.mkdir(exist_ok=True)

# Create config.toml with dark theme
config_content = '''[theme]
primaryColor = "#FF4B4B"
backgroundColor = "#0E1117"
secondaryBackgroundColor = "#262730"
textColor = "#FAFAFA"
font = "sans serif"

[server]
headless = true
'''

config_file = config_dir / 'config.toml'
config_file.write_text(config_content, encoding='utf-8')

print('✅ Created dark theme configuration')
print(f'📁 Location: {config_file.absolute()}')
print()
print('🎨 Theme colors:')
print('  • Background: Dark gray (#0E1117)')
print('  • Secondary: Darker gray (#262730)')
print('  • Primary: Red (#FF4B4B)')
print('  • Text: White (#FAFAFA)')
print()
print('🔄 Restart your app to see dark mode!')
