use std::path::Path;
use std::str::FromStr;

// ---------------------------------------------------------------------------
// Tag category enum
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum TagCategory {
    Even,
    Odd,
    Y,
    Rpm,
    Dpm,
    LigTag,
    Spacer,
}

impl FromStr for TagCategory {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_uppercase().trim() {
            "EVEN"   => Ok(TagCategory::Even),
            "ODD"    => Ok(TagCategory::Odd),
            "Y"      => Ok(TagCategory::Y),
            "RPM"    => Ok(TagCategory::Rpm),
            "DPM"    => Ok(TagCategory::Dpm),
            "LIGTAG" => Ok(TagCategory::LigTag),
            "SPACER" => Ok(TagCategory::Spacer),
            other    => Err(format!("Unknown tag category: '{}'", other)),
        }
    }
}

impl std::fmt::Display for TagCategory {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = match self {
            TagCategory::Even   => "EVEN",
            TagCategory::Odd    => "ODD",
            TagCategory::Y      => "Y",
            TagCategory::Rpm    => "RPM",
            TagCategory::Dpm    => "DPM",
            TagCategory::LigTag => "LIGTAG",
            TagCategory::Spacer => "SPACER",
        };
        write!(f, "{}", s)
    }
}

// ---------------------------------------------------------------------------
// Tag definition (one entry from the config file)
// ---------------------------------------------------------------------------

#[derive(Debug, Clone)]
pub struct TagDef {
    pub category:   TagCategory,
    pub name:       String,
    pub seq:        Vec<u8>,
    pub mismatches: u8,
}

// ---------------------------------------------------------------------------
// Config
// ---------------------------------------------------------------------------

const DEFAULT_SPACER: usize = 6;
const DEFAULT_LAXITY: usize = 6;

#[derive(Debug, Clone)]
pub struct Config {
    pub layout1:    Vec<TagCategory>,
    pub layout2:    Vec<TagCategory>,
    pub spacer_len: usize,
    pub laxity:     usize,
    pub tag_defs:   Vec<TagDef>,
}

impl Config {
    pub fn from_file(path: &Path) -> Result<Self, Box<dyn std::error::Error>> {
        let content = std::fs::read_to_string(path)?;
        Self::from_str_content(&content)
    }

    /// Parse config from a string — used directly in unit tests.
    pub fn from_str_content(content: &str) -> Result<Self, Box<dyn std::error::Error>> {
        let mut layout1: Vec<TagCategory> = Vec::new();
        let mut layout2: Vec<TagCategory> = Vec::new();
        let mut spacer_len: Option<usize> = None;
        let mut laxity: Option<usize> = None;
        let mut tag_defs: Vec<TagDef> = Vec::new();

        for (line_no, raw_line) in content.lines().enumerate() {
            let line = raw_line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            if line.starts_with("READ1 = ") {
                if !layout1.is_empty() {
                    log::warn!("LAYOUT found twice in config (line {}). Overwriting.", line_no + 1);
                    layout1.clear();
                }
                parse_layout(&mut layout1, line)?;
            } else if line.starts_with("READ2 = ") {
                if !layout2.is_empty() {
                    log::warn!("LAYOUT found twice in config (line {}). Overwriting.", line_no + 1);
                    layout2.clear();
                }
                parse_layout(&mut layout2, line)?;
            } else if line.starts_with("SPACER = ") {
                if spacer_len.is_some() {
                    log::warn!("SPACER found twice in config. Overwriting old value.");
                }
                let val = line["SPACER = ".len()..].trim().parse::<usize>()
                    .map_err(|e| format!("Bad SPACER value: {}", e))?;
                spacer_len = Some(val);
            } else if line.starts_with("LAXITY = ") {
                if laxity.is_some() {
                    log::warn!("LAXITY found twice in config. Overwriting old value.");
                }
                let val = line["LAXITY = ".len()..].trim().parse::<usize>()
                    .map_err(|e| format!("Bad LAXITY value: {}", e))?;
                laxity = Some(val);
            } else {
                // Tag line: CATEGORY\tNAME\tSEQUENCE\tNUM_MISMATCHES
                match parse_tag_line(line) {
                    Ok(def) => tag_defs.push(def),
                    Err(e) => {
                        return Err(format!("Error on config line {}: {}", line_no + 1, e).into())
                    }
                }
            }
        }

        if layout1.is_empty() && layout2.is_empty() {
            return Err("No tag layouts found in configuration file.".into());
        }
        if tag_defs.is_empty() {
            return Err("No tags found in configuration file.".into());
        }

        Ok(Config {
            layout1,
            layout2,
            spacer_len: spacer_len.unwrap_or(DEFAULT_SPACER),
            laxity:     laxity.unwrap_or(DEFAULT_LAXITY),
            tag_defs,
        })
    }
}

fn parse_layout(
    layout: &mut Vec<TagCategory>,
    line: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    // "READ1 = DPM"  or  "READ2 = Y|SPACER|ODD|SPACER|EVEN"
    let rhs = line
        .splitn(2, "= ")
        .nth(1)
        .ok_or("Malformed READ line (missing '= ')")?;
    for part in rhs.split('|') {
        layout.push(part.trim().parse::<TagCategory>()
            .map_err(|e| format!("In layout '{}': {}", rhs, e))?);
    }
    Ok(())
}

fn parse_tag_line(line: &str) -> Result<TagDef, Box<dyn std::error::Error>> {
    let fields: Vec<&str> = line.splitn(4, '\t').collect();
    if fields.len() != 4 {
        return Err(format!(
            "Expected 4 tab-separated fields, got {}: {:?}",
            fields.len(),
            line
        ).into());
    }
    let category:   TagCategory = fields[0].parse()?;
    let name:       String      = fields[1].trim().to_string();
    let seq:        Vec<u8>     = fields[2].trim().as_bytes().to_vec();
    let mismatches: u8          = fields[3].trim().parse()
        .map_err(|e| format!("Bad mismatch count '{}': {}", fields[3].trim(), e))?;

    Ok(TagDef { category, name, seq, mismatches })
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const MINI_CONFIG: &str = "\
READ1 = DPM
READ2 = Y|SPACER|ODD
SPACER = 3
LAXITY = 2
DPM\tDPM_A\tAAAACCCC\t0
DPM\tDPM_B\tGGGGTTTT\t0
Y\tY_1\tTTTTAAAA\t0
Y\tY_2\tCCCCGGGG\t0
ODD\tOdd_1\tATCGATCG\t1
ODD\tOdd_2\tGCTAGCTA\t1
";

    #[test]
    fn test_parse_layouts() {
        let cfg = Config::from_str_content(MINI_CONFIG).unwrap();
        assert_eq!(cfg.layout1, vec![TagCategory::Dpm]);
        assert_eq!(
            cfg.layout2,
            vec![TagCategory::Y, TagCategory::Spacer, TagCategory::Odd]
        );
    }

    #[test]
    fn test_parse_spacer_and_laxity() {
        let cfg = Config::from_str_content(MINI_CONFIG).unwrap();
        assert_eq!(cfg.spacer_len, 3);
        assert_eq!(cfg.laxity, 2);
    }

    #[test]
    fn test_spacer_default() {
        let content = "READ1 = DPM\nDPM\tDPM_A\tAAAA\t0\n";
        let cfg = Config::from_str_content(content).unwrap();
        assert_eq!(cfg.spacer_len, DEFAULT_SPACER);
    }

    #[test]
    fn test_laxity_default() {
        let content = "READ1 = DPM\nDPM\tDPM_A\tAAAA\t0\n";
        let cfg = Config::from_str_content(content).unwrap();
        assert_eq!(cfg.laxity, DEFAULT_LAXITY);
    }

    #[test]
    fn test_tag_count() {
        let cfg = Config::from_str_content(MINI_CONFIG).unwrap();
        assert_eq!(cfg.tag_defs.len(), 6);
    }

    #[test]
    fn test_tag_fields() {
        let cfg = Config::from_str_content(MINI_CONFIG).unwrap();
        let dpm_a = cfg.tag_defs.iter().find(|t| t.name == "DPM_A").unwrap();
        assert_eq!(dpm_a.category, TagCategory::Dpm);
        assert_eq!(dpm_a.seq, b"AAAACCCC");
        assert_eq!(dpm_a.mismatches, 0);
    }

    #[test]
    fn test_comments_and_blanks_ignored() {
        let content = "# comment\n\nREAD1 = DPM\n# another comment\nDPM\tA\tATCG\t0\n";
        let cfg = Config::from_str_content(content).unwrap();
        assert_eq!(cfg.tag_defs.len(), 1);
    }

    #[test]
    fn test_no_layouts_error() {
        let content = "DPM\tA\tATCG\t0\n";
        assert!(Config::from_str_content(content).is_err());
    }

    #[test]
    fn test_no_tags_error() {
        let content = "READ1 = DPM\n";
        assert!(Config::from_str_content(content).is_err());
    }

    #[test]
    fn test_unknown_category_error() {
        let content = "READ1 = UNKNOWN\nDPM\tA\tATCG\t0\n";
        assert!(Config::from_str_content(content).is_err());
    }

    #[test]
    fn test_category_case_insensitive() {
        let content = "READ1 = dpm\ndpm\tA\tATCG\t0\n";
        let cfg = Config::from_str_content(content).unwrap();
        assert_eq!(cfg.layout1[0], TagCategory::Dpm);
    }
}
