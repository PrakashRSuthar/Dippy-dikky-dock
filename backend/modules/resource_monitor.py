# backend/modules/standalone_resource_monitor.py
import psutil
import time
from datetime import datetime

class StandaloneResourceMonitor:
    """
    Simple resource monitor that runs independently.
    Start it before docking, stop it manually after docking is done.
    """
    
    def __init__(self, interval=5):
        self.interval = interval
        self.data = []
        self.start_time = None
    
    def monitor(self):
        """Start monitoring - runs until Ctrl+C"""
        self.start_time = datetime.now()
        print(f"🔍 Resource monitoring started at {self.start_time.strftime('%H:%M:%S')}")
        print(f"📊 Sampling every {self.interval} seconds")
        print("🚀 Start your docking pipeline now")
        print("⏹️  Press Ctrl+C to stop monitoring and get report")
        print("-" * 50)
        
        try:
            sample_count = 0
            while True:
                # Get current resource usage
                cpu = psutil.cpu_percent(interval=1)
                memory = psutil.virtual_memory()
                disk = psutil.disk_usage('.')
                
                # Store snapshot
                snapshot = {
                    'timestamp': datetime.now(),
                    'cpu_percent': cpu,
                    'memory_percent': memory.percent,
                    'memory_used_gb': (memory.total - memory.available) / (1024**3),
                    'memory_available_gb': memory.available / (1024**3),
                    'disk_percent': (disk.used / disk.total) * 100,
                    'disk_free_gb': disk.free / (1024**3)
                }
                
                self.data.append(snapshot)
                sample_count += 1
                
                # Show live status every 10 samples
                if sample_count % 10 == 0:
                    elapsed = (datetime.now() - self.start_time).total_seconds() / 60
                    print(f"[{elapsed:05.1f}min] CPU: {cpu:5.1f}% | Memory: {snapshot['memory_used_gb']:4.1f}GB | Free Disk: {snapshot['disk_free_gb']:5.1f}GB")
                
                time.sleep(self.interval - 1)  # -1 because cpu_percent takes 1 second
                
        except KeyboardInterrupt:
            self.generate_report()
    
    def generate_report(self):
        """Generate and display resource usage report"""
        if not self.data:
            print("❌ No data collected")
            return
        
        end_time = datetime.now()
        duration = (end_time - self.start_time).total_seconds() / 60
        
        # Calculate statistics
        cpu_values = [d['cpu_percent'] for d in self.data]
        memory_used = [d['memory_used_gb'] for d in self.data]
        memory_available = [d['memory_available_gb'] for d in self.data]
        disk_free = [d['disk_free_gb'] for d in self.data]
        
        print(f"\n" + "="*60)
        print(f"📊 RESOURCE USAGE REPORT")
        print(f"="*60)
        
        print(f"⏱️  Monitoring Duration: {duration:.1f} minutes")
        print(f"📈 Total Samples: {len(self.data)}")
        print(f"🕐 Start Time: {self.start_time.strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"🕐 End Time: {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
        
        print(f"\n💻 CPU Usage:")
        print(f"   • Average: {sum(cpu_values)/len(cpu_values):.1f}%")
        print(f"   • Minimum: {min(cpu_values):.1f}%")
        print(f"   • Maximum: {max(cpu_values):.1f}%")
        
        print(f"\n🧠 Memory Usage:")
        print(f"   • Average Used: {sum(memory_used)/len(memory_used):.1f} GB")
        print(f"   • Average Available: {sum(memory_available)/len(memory_available):.1f} GB")
        print(f"   • Peak Used: {max(memory_used):.1f} GB")
        print(f"   • Minimum Available: {min(memory_available):.1f} GB")
        
        print(f"\n💾 Disk Usage:")
        print(f"   • Average Free Space: {sum(disk_free)/len(disk_free):.1f} GB")
        print(f"   • Minimum Free Space: {min(disk_free):.1f} GB")
        
        # Performance analysis
        avg_cpu = sum(cpu_values) / len(cpu_values)
        avg_memory = sum(memory_used) / len(memory_used)
        
        print(f"\n🔍 Performance Analysis:")
        if avg_cpu < 50:
            print(f"   • CPU: Low usage ({avg_cpu:.1f}%) - System can handle more parallel jobs")
        elif avg_cpu < 80:
            print(f"   • CPU: Moderate usage ({avg_cpu:.1f}%) - Good utilization")
        else:
            print(f"   • CPU: High usage ({avg_cpu:.1f}%) - Consider reducing parallel jobs")
        
        if avg_memory < 4:
            print(f"   • Memory: Low usage ({avg_memory:.1f}GB) - Can handle more ligands")
        elif avg_memory < 8:
            print(f"   • Memory: Moderate usage ({avg_memory:.1f}GB) - Good utilization")
        else:
            print(f"   • Memory: High usage ({avg_memory:.1f}GB) - Consider smaller batches")
        
        # Recommendations for system check
        system_cores = psutil.cpu_count()
        if avg_cpu < 60:
            recommended_parallel = min(system_cores, int(system_cores * 0.8))
        else:
            recommended_parallel = min(system_cores // 2, 4)
        
        print(f"\n💡 Recommendations for System Check Module:")
        print(f"   • Suggested max parallel ligands: {recommended_parallel}")
        print(f"   • Estimated memory per ligand: {avg_memory:.1f} GB (if running parallel)")
        print(f"   • System utilization: {'Good' if 50 <= avg_cpu <= 80 else 'Adjust needed'}")
        
        print(f"="*60)


if __name__ == "__main__":
    print("🚀 Standalone Resource Monitor")
    print("This will monitor your system independently")
    print("Start this, then run docking in another terminal")
    print()
    
    # Create and run monitor
    monitor = StandaloneResourceMonitor(interval=5)
    monitor.monitor()
