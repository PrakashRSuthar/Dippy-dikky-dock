# backend/modules/simple_system_check.py
import multiprocessing
import psutil
import shutil

class SimpleSystemChecker:
    """Simple system check that gives max safe ligand limit upfront"""
    
    def __init__(self):
        # Resource requirements per ligand
        self.memory_per_ligand = 0.5  # GB
        self.disk_per_ligand = 0.1    # GB
    
    def check_system_limits(self):
        """Check system and return max safe parallel ligands"""
        
        # Get system specs
        cpu_cores = multiprocessing.cpu_count()
        memory_free_gb = psutil.virtual_memory().available / (1024**3)
        disk_free_gb = shutil.disk_usage('.').free / (1024**3)
        
        # Calculate conservative limits
        max_by_cpu = max(1, cpu_cores - 1)  # Keep 1 core free
        max_by_memory = max(1, int(memory_free_gb / self.memory_per_ligand))
        max_by_disk = max(1, int(disk_free_gb / self.disk_per_ligand))
        
        # Take most restrictive limit
        max_safe_parallel = min(max_by_cpu, max_by_memory, max_by_disk)
        
        return {
            'cpu_cores': cpu_cores,
            'memory_gb': round(memory_free_gb, 1),
            'disk_gb': round(disk_free_gb, 1),
            'max_parallel': max_safe_parallel
        }
    
    def display_system_check(self):
        """Display system specs and limits to user"""
        
        specs = self.check_system_limits()
        
        print("🖥️  SYSTEM CHECK")
        print("="*40)
        print(f"CPU Cores: {specs['cpu_cores']}")
        print(f"Available Memory: {specs['memory_gb']} GB")
        print(f"Free Disk: {specs['disk_gb']} GB")
        print("="*40)
        print(f"💊 MAX SAFE PARALLEL LIGANDS: {specs['max_parallel']}")
        
        if specs['max_parallel'] > 1:
            print("✅ Your PC can handle parallel processing")
        else:
            print("⚠️  Your PC should process one ligand at a time")
        
        print("="*40)
        
        return specs['max_parallel']


def get_system_max_limit():
    """Quick function to get max safe parallel limit"""
    
    print("🔍 Checking your system capabilities...")
    
    checker = SimpleSystemChecker()
    max_limit = checker.display_system_check()
    
    return max_limit


# Usage in your pipeline
def batch_setup_with_limit_check():
    """Setup batch processing with automatic limit check"""
    
    # Check system first
    max_limit = get_system_max_limit()
    
    # Then ask user for ligands with limit shown
    while True:
        try:
            user_input = input(f"\nHow many ligands do you want to test? (Max recommended: {max_limit}): ")
            num_ligands = int(user_input)
            
            if num_ligands <= 0:
                print("Please enter a positive number")
                continue
            elif num_ligands > max_limit:
                print(f"⚠️  Warning: {num_ligands} exceeds safe limit of {max_limit}")
                confirm = input("Continue anyway? This may slow down your PC [y/n]: ")
                if confirm.lower() != 'y':
                    continue
            
            # Determine processing mode
            if max_limit > 1 and num_ligands > 1:
                batch_size = min(max_limit, num_ligands)
                mode = "parallel"
                print(f"✅ Will process {batch_size} ligands at once")
            else:
                batch_size = 1
                mode = "sequential"
                print(f"✅ Will process 1 ligand at a time")
            
            return {
                'num_ligands': num_ligands,
                'batch_size': batch_size,
                'mode': mode,
                'max_limit': max_limit
            }
            
        except ValueError:
            print("Please enter a valid number")
        except KeyboardInterrupt:
            print("\n❌ Cancelled")
            return None


if __name__ == "__main__":
    config = batch_setup_with_limit_check()
    
    if config:
        print(f"\n🚀 Configuration:")
        print(f"Ligands to test: {config['num_ligands']}")
        print(f"Batch size: {config['batch_size']}")
        print(f"Mode: {config['mode']}")
    else:
        print("Setup cancelled")
