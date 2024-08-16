using System;
using System.Collections.Generic;
using System.ComponentModel.DataAnnotations;
using System.ComponentModel.DataAnnotations.Schema;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Infrastructure.Models
{
    public class Block
    {
        [Key]
        public int Id { get; set; }

        [Key]
        public Rental? Rental { get; set; }

        public string? BlockedStart { get; set; }
        public string? BlockedEnd { get; set; }
    }
}
